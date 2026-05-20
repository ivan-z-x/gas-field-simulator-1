from __future__ import annotations
import math
import numpy as np
from src.fluid import Fluid
from src.pipe import Pipe

class Well:
    def __init__(self, fluid: Fluid, k: float, h: float, re: float, rw: float, pipe: Pipe | None = None,
                 skin: float = 0.0):
        self.fluid = fluid
        self.k = k          # мД
        self.h = h          # м (эффективная мощность пласта)
        self.re = re        # м (радиус дренирования)
        self.rw = rw        # м (радиус скважины)
        self.skin = skin    # скин-фактор
        self.pipe = pipe

    def productivity_index(self, P_res: float) -> float:
        """Индекс продуктивности (линейный закон Дарси для газа с приближением)."""
        mu = self.fluid.mu(P_res)   # сП
        # Коэффициент пропорциональности подобран, чтобы дебиты были порядка тысяч м³/сут
        C = (self.k * self.h) / (mu * math.log(self.re / self.rw)) * 0.01
        return max(0.0, C)

    def q(self, P_res: float, P_bhp: float) -> float:
        """Дебит скважины при заданных пластовом и забойном давлениях."""
        if P_bhp >= P_res:
            return 0.0
        return self.productivity_index(P_res) * (P_res - P_bhp)

    def ipr_curve(self, P_res: float, n: int = 50):
        """Кривая IPR (зависимость дебита от забойного давления)."""
        xs = np.linspace(0.0, P_res, n)
        ys = [self.q(P_res, p) for p in xs]
        return xs, np.array(ys)

    # ---------- Добавленная функция VLP для гидравлики НКТ ----------
    def vlp(self, P_wh: float, q_std: float) -> float:
        """
        Рассчитывает забойное давление P_bhp по устьевому давлению P_wh и дебиту q_std.
        Использует self.pipe (НКТ).
        """
        if q_std <= 0:
            # Гидростатический столб газа (приближённо)
            rho_avg = self.fluid.ro(P_wh)
            dP_hyd = rho_avg * 9.81 * self.pipe.vertical_depth / 101325.0
            return P_wh + dP_hyd
        P_bhp = P_wh + 5.0
        for _ in range(10):
            state = self.pipe.dp(P_bhp, q_std, name="vlp")
            if state.P_out <= 0:
                return P_wh
            if abs(state.P_out - P_wh) < 1e-3:
                return state.P_in
            P_bhp = P_bhp + (P_wh - state.P_out) * 0.5
            if P_bhp < P_wh:
                P_bhp = P_wh + 1.0
        return P_bhp
    # -------------------------------------------------------------