from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from src.state import NodeState
from src.reservoir import Reservoir
from src.pipe import Pipe
from src.well import Well
from src.compressor import DCS


class FieldSimulator:
    """Симулятор газового куста: три скважины на общем пласте."""

    def __init__(self, reservoir: Reservoir, wells: list[Well],
                 shlyf: Pipe, dcs: DCS):
        if len(wells) != 3:
            raise ValueError("Требуется ровно 3 скважины")
        self.reservoir = reservoir
        self.wells     = wells
        self.shlyf     = shlyf
        self.dcs       = dcs

    def solve(self, P_res: float) -> dict[str, NodeState]:
        """Нахождение рабочей точки при заданном P_res [атм]."""
        P_res  = float(P_res)
        P_dcs  = self.dcs.P_in()

        def residuals(x):
            q1, q2, q3, P_man = x[0], x[1], x[2], x[3]
            F = []
            for qi, well in zip([q1, q2, q3], self.wells):
                qi = max(0.0, float(qi))
                dp_tube = well.pipe.dp(max(P_dcs, float(P_man)), qi).dP if qi > 1e-6 else 0.0
                P_bhp   = float(P_man) + dp_tube
                F.append(qi - well.q(P_res, P_bhp))
            q_total = max(0.0, q1) + max(0.0, q2) + max(0.0, q3) + self.dcs.q_ext
            dp_shlyf = self.shlyf.dp(P_dcs, max(1.0, q_total)).dP
            F.append(float(P_man) - (P_dcs + dp_shlyf))
            return F

        # Начальное приближение
        C1 = self.wells[0].productivity_index(P_res)
        q0 = max(1.0, C1 * (P_res - P_dcs) * 0.5)
        x0 = np.array([q0, q0, q0, P_dcs + 0.1])
        sol = fsolve(residuals, x0, full_output=True)
        x   = sol[0]

        q1 = max(0.0, float(x[0]))
        q2 = max(0.0, float(x[1]))
        q3 = max(0.0, float(x[2]))
        P_man = max(P_dcs, float(x[3]))

        states: dict[str, NodeState] = {}
        for idx, (qi, well) in enumerate(zip([q1, q2, q3], self.wells), 1):
            dp = well.pipe.dp(P_man, qi).dP if qi > 1e-6 else 0.0
            st = well.pipe.dp(P_man, qi) if qi > 1e-6 else None
            states[f"well_{idx}"] = NodeState(
                name=f"well_{idx}", P_in=P_man+dp, P_out=P_man,
                dP=dp, q_std=qi,
                q_res=st.q_res if st else 0.0,
                v=st.v if st else 0.0,
                rho=st.rho if st else None,
            )

        q_total = q1 + q2 + q3 + self.dcs.q_ext
        st_sh   = self.shlyf.dp(P_dcs, max(1.0, q_total))
        states["shlyf"] = NodeState(
            name="shlyf", P_in=P_man, P_out=P_dcs,
            dP=st_sh.dP, q_std=q_total,
            q_res=st_sh.q_res, v=st_sh.v, rho=st_sh.rho,
        )
        states["dcs"] = NodeState(
            name="dcs", P_in=P_dcs, P_out=self.dcs.P_line,
            dP=self.dcs.P_line - P_dcs, q_std=q_total,
        )
        return states

    def run(self, N_days: int, dt: float = 1.0) -> pd.DataFrame:
        """Динамическая симуляция на N_days шагов с шагом dt [сут]."""
        rows = []
        for t in range(N_days):
            P_res  = float(self.reservoir.resprops.P)
            states = self.solve(P_res)
            q1 = float(states["well_1"].q_std)
            q2 = float(states["well_2"].q_std)
            q3 = float(states["well_3"].q_std)
            q_wells = q1 + q2 + q3
            P_new   = self.reservoir.p2(q_wells, dt)
            self.reservoir.resprops.P = P_new
            rows.append({
                "t, сут":             t,
                "P_res, атм":         P_res,
                "P_man, атм":         float(states["shlyf"].P_in),
                "q1, ст.м³/сут":      q1,
                "q2, ст.м³/сут":      q2,
                "q3, ст.м³/сут":      q3,
                "q_total, ст.м³/сут": q_wells,
                "Gp, тыс.ст.м³":      0.0,
            })
        df = pd.DataFrame(rows)
        df["Gp, тыс.ст.м³"] = (df["q_total, ст.м³/сут"] * dt).cumsum() / 1000.0
        return df
