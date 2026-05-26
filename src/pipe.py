from __future__ import annotations
import math
from src.fluid import Fluid
from src.state import NodeState


class Pipe:
    """
    Гидравлическая модель трубопровода (НКТ или шлейф).
    Уравнение Дарси–Вейсбаха + Колбрук–Уайт.
    """

    def __init__(self, L: float, D: float, roughness: float,
                 fluid: Fluid, vertical_depth: float = 0.0):
        """
        L             — длина трубы [м]
        D             — внутренний диаметр [м]
        roughness (δ) — абсолютная шероховатость [м]
        vertical_depth (H) — вертикальная глубина [м] (0 для шлейфа)
        """
        self.L = L
        self.D = D
        self.roughness = roughness
        self.fluid = fluid
        self.vertical_depth = vertical_depth

    def _velocity(self, P: float, q_std: float) -> float:
        """Средняя скорость газа в трубе [м/с]."""
        bg = self.fluid.bg(P)
        area = math.pi * self.D**2 / 4.0
        # q_std [ст.м³/сут] → [м³/с]: делим на 86400
        return q_std * bg / (area * 86400.0)

    def _reynolds(self, P: float, q_std: float) -> float:
        """Число Рейнольдса."""
        rho   = self.fluid.ro(P)
        v     = self._velocity(P, q_std)
        mu_pa = self.fluid.mu(P) / 1000.0   # сП → Па·с
        return rho * v * self.D / max(mu_pa, 1e-12)

    def _friction_factor(self, Re: float) -> float:
        """Коэффициент гидравлического сопротивления λ (Колбрук–Уайт)."""
        if Re <= 0:
            return 0.0
        if Re < 2300:
            return 64.0 / Re   # Пуазейль (ламинарный режим)
        # Итерации Колбрука–Уайта: λ^(n+1) = [-2 lg(δ/3.7D + 2.51/Re√λ)]^(-2)
        lam = 0.02
        for _ in range(100):
            arg = self.roughness / (3.7 * self.D) + 2.51 / (Re * math.sqrt(lam))
            lam_new = 1.0 / (2.0 * math.log10(arg))**2
            if abs(lam_new - lam) < 1e-6:
                return lam_new
            lam = lam_new
        return lam

    def dp(self, P: float, q: float, name: str = "pipe") -> NodeState:
        """
        Перепад давления Дарси–Вейсбаха.
        P [атм] — давление на входе, q [ст.м³/сут].
        Возвращает NodeState: P_in=входное, P_out=выходное давление.
        """
        P = float(P)
        q = float(q)
        rho = self.fluid.ro(P)
        bg  = self.fluid.bg(P)
        v   = self._velocity(P, q)
        Re  = self._reynolds(P, q)
        lam = self._friction_factor(Re)

        # ΔP = λ·L/D·ρv²/2 + ρgH  [Па]
        dp_pa  = lam * (self.L / self.D) * rho * v**2 / 2.0
        dp_pa += rho * 9.81 * self.vertical_depth
        dp_atm = dp_pa / 101325.0

        q_res = q * bg   # [м³/сут]
        return NodeState(
            name=name,
            P_in=P,
            P_out=max(0.0, P - dp_atm),
            dP=dp_atm,
            q_std=q,
            q_res=q_res,
            v=v,
            rho=rho,
        )
