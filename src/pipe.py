from __future__ import annotations
import math
from src.fluid import Fluid
from src.state import NodeState

class Pipe:
    def __init__(self, L: float, D: float, roughness: float, fluid: Fluid, vertical_depth: float = 0.0):
        self.L = L
        self.D = D
        self.roughness = roughness
        self.fluid = fluid
        self.vertical_depth = vertical_depth

    def q_res(self, P: float, q_std: float) -> float:
        """Перевод дебита в пластовые условия через Bg."""
        return q_std * self.fluid.bg(P)

    def velocity(self, P: float, q_std: float) -> float:
        """Скорость потока в трубе [м/с]."""
        q_res_day = self.q_res(P, q_std)
        area = math.pi * self.D**2 / 4.0
        return q_res_day / 86400.0 / area

    def reynolds(self, P: float, q_std: float) -> float:
        """Число Рейнольдса."""
        rho = self.fluid.ro(P)
        v = self.velocity(P, q_std)
        mu_pas = self.fluid.mu(P) / 1000.0  # Па·с
        return rho * v * self.D / max(mu_pas, 1e-12)

    def friction_factor(self, Re: float) -> float:
        """Коэффициент трения (ламинарный или итерация Кольбрука-Уайта)."""
        if Re <= 0:
            return 0.0
        if Re < 2300:
            return 64.0 / Re
        lam = 0.02
        for _ in range(100):
            rhs = -2.0 * math.log10(self.roughness / (3.7 * self.D) + 2.51 / (Re * math.sqrt(lam)))
            lam_new = 1.0 / (rhs**2)
            if abs(lam_new - lam) < 1e-6:
                return lam_new
            lam = lam_new
        return lam

    def dp(self, P: float, q: float, name: str = "pipe") -> NodeState:
        """Расчёт перепада давления на трубе (включая гидростатику)."""
        rho = self.fluid.ro(P)
        q_res = self.q_res(P, q)
        v = self.velocity(P, q)
        Re = self.reynolds(P, q)
        lam = self.friction_factor(Re)
        # Потери на трение + гидростатика
        dp_pa = lam * (self.L / self.D) * rho * v**2 / 2.0 + rho * 9.81 * self.vertical_depth
        dp_atm = dp_pa / 101325.0
        return NodeState(
            name=name,
            P_in=float(P),
            P_out=max(0.0, float(P) - dp_atm),
            dP=dp_atm,
            q_std=float(q),
            q_res=q_res,
            v=v,
            rho=rho,
        )