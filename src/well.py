from __future__ import annotations
import math
from src.fluid import Fluid
from src.pipe import Pipe

# Коэффициент перевода единиц Дарси (мД, м, атм, м³/сут)
BETA = 0.00852702   # ст.м³/(сут·мД·м·атм)


class Well:
    """
    Модель скважины: приток газа из пласта (закон Дарси) + НКТ (Pipe).
    """

    def __init__(self, fluid: Fluid, k: float, h: float,
                 re: float, rw: float, pipe: Pipe | None = None,
                 skin: float = 0.0):
        """
        k    — проницаемость [мД]
        h    — эффективная мощность [м]
        re   — радиус контура питания [м]
        rw   — радиус скважины [м]
        pipe — НКТ (объект Pipe)
        skin — скин-фактор
        """
        self.fluid = fluid
        self.k = k
        self.h = h
        self.re = re
        self.rw = rw
        self.skin = skin
        self.pipe = pipe

    def productivity_index(self, P_res: float) -> float:
        """
        Коэффициент продуктивности C [ст.м³/(сут·атм)].
        C = β·k·h / (μ·ln(re/rw) + skin)
        """
        mu = self.fluid.mu(P_res)   # сП
        ln_term = math.log(self.re / self.rw) + self.skin
        return BETA * self.k * self.h / (mu * ln_term)

    def q(self, P_res: float, P_bhp: float) -> float:
        """
        Дебит скважины по закону Дарси [ст.м³/сут].
        q = C · (P_res - P_bhp)
        """
        if P_bhp >= P_res:
            return 0.0
        return self.productivity_index(P_res) * (P_res - P_bhp)

    def ipr_curve(self, P_res: float, n: int = 50):
        """Кривая IPR: списки (P_bhp, q)."""
        step = P_res / (n - 1)
        P_bhp_vals = [i * step for i in range(n)]
        q_vals     = [self.q(P_res, p) for p in P_bhp_vals]
        return P_bhp_vals, q_vals

    def vlp_curve(self, P_man: float, n: int = 50):
        """
        Кривая VLP при фиксированном манифолдном давлении.
        Возвращает (q_vals, P_bhp_vals).
        P_bhp = P_man + ΔP_НКТ(q)
        """
        q_max = self.q(self.fluid.T, 0.0) if False else 5000.0
        step  = q_max / (n - 1)
        q_vals     = [i * step for i in range(n)]
        P_bhp_vals = []
        for qi in q_vals:
            if qi <= 0:
                P_bhp_vals.append(P_man)
            else:
                state = self.pipe.dp(P_man + self.pipe.dp(P_man, qi).dP, qi)
                # P_bhp = P_man + ΔP_НКТ
                st = self.pipe.dp(P_man, qi)
                P_bhp_vals.append(P_man + st.dP)
        return q_vals, P_bhp_vals
