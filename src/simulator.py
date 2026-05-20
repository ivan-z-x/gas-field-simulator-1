from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import root
from src.state import NodeState
from src.reservoir import Reservoir
from src.pipe import Pipe
from src.well import Well
from src.compressor import DCS

class FieldSimulator:
    def __init__(self, reservoir: Reservoir, wells: list[Well], shlyf: Pipe, dcs: DCS):
        if len(wells) != 3:
            raise ValueError("Требуется ровно 3 скважины")
        self.reservoir = reservoir
        self.wells = wells
        self.shlyf = shlyf
        self.dcs = dcs

    def solve(self, P_res: float) -> dict[str, NodeState]:
        """
        Решение системы нелинейных уравнений для рабочей точки.
        """
        def residuals(x):
            q1, q2, q3, P_man = x
            qs = [max(0.0, float(q1)), max(0.0, float(q2)), max(0.0, float(q3))]
            eqs = []
            for qi, well in zip(qs, self.wells):
                if qi <= 0:
                    eqs.append(qi)
                    continue
                tube = well.pipe.dp(P_man + 1e-9, qi, name="tube")
                P_bhp = P_man + tube.dP
                q_ipr = well.q(P_res, P_bhp)
                eqs.append(qi - q_ipr)
            q_total = sum(qs) + self.dcs.q_ext
            shlyf = self.shlyf.dp(self.dcs.P_in() + 1e-9, q_total, name="shlyf")
            eqs.append(P_man - (self.dcs.P_in() + shlyf.dP))
            return np.array(eqs, dtype=float)

        x0 = np.array([5000.0, 5000.0, 5000.0, self.dcs.P_in() + 5.0], dtype=float)
        sol = root(residuals, x0, method='hybr', options={'maxfev': 2000})
        x = sol.x
        q1 = max(0.0, float(x[0]))
        q2 = max(0.0, float(x[1]))
        q3 = max(0.0, float(x[2]))
        P_man = max(0.0, float(x[3]))

        states = {}
        q_list = [q1, q2, q3]
        for idx, (qi, well) in enumerate(zip(q_list, self.wells), 1):
            tube = well.pipe.dp(P_man + 1e-9, qi, name=f"well_{idx}")
            P_bhp = P_man + tube.dP
            # Приведение qi к скалярному float
            if isinstance(qi, (np.ndarray, list, tuple)):
                qi = float(qi[0]) if len(qi) > 0 else 0.0
            else:
                qi = float(qi)
            states[f"well_{idx}"] = NodeState(
                name=f"well_{idx}",
                P_in=float(P_bhp),
                P_out=float(P_man),
                dP=float(tube.dP),
                q_std=qi,
                q_res=tube.q_res,
                v=tube.v,
                rho=tube.rho,
            )
        q_total = q1 + q2 + q3 + self.dcs.q_ext
        shlyf = self.shlyf.dp(self.dcs.P_in() + 1e-9, q_total, name="shlyf")
        states["shlyf"] = NodeState(
            name="shlyf",
            P_in=float(P_man),
            P_out=float(self.dcs.P_in()),
            dP=float(shlyf.dP),
            q_std=float(q_total),
            q_res=shlyf.q_res,
            v=shlyf.v,
            rho=shlyf.rho,
        )
        states["dcs"] = NodeState(
            name="dcs",
            P_in=float(self.dcs.P_in()),
            P_out=float(self.dcs.P_line),
            dP=float(self.dcs.P_line - self.dcs.P_in()),
            q_std=float(q_total),
            q_res=None,
            v=None,
            rho=None,
        )
        return states

    def run(self, N_days: int, dt: float = 1.0) -> pd.DataFrame:
        """Запуск динамической симуляции на N_days дней с шагом dt."""
        rows = []
        self.reservoir.Gp = 0.0
        for t in range(N_days):
            P_res = self.reservoir.resprops.P
            states = self.solve(P_res)
            q1 = states["well_1"].q_std
            q2 = states["well_2"].q_std
            q3 = states["well_3"].q_std
            # Приведение к скаляру
            def to_scalar(v):
                if isinstance(v, (np.ndarray, list, tuple)):
                    return float(v[0]) if len(v) > 0 else 0.0
                return float(v)
            q1 = to_scalar(q1)
            q2 = to_scalar(q2)
            q3 = to_scalar(q3)
            q_total = q1 + q2 + q3
            self.reservoir.Gp += q_total * dt
            Gp_thousand = self.reservoir.Gp / 1000.0
            rows.append({
                "t, day": t,
                "P_res, atm": float(P_res),
                "P_man, atm": float(states["shlyf"].P_in),
                "q1, std m3/day": q1,
                "q2, std m3/day": q2,
                "q3, std m3/day": q3,
                "q_total, std m3/day": q_total,
                "Gp, thousand std m3": Gp_thousand,
            })
            new_P = self.reservoir.p2(q_total, dt)
            self.reservoir.resprops.P = new_P
            if (t+1) % 10 == 0 or t == N_days-1:
                print(f"Шаг {t+1}/{N_days}: P_пл={new_P:.2f} атм, q_сум={q_total:.1f} ст.м³/сут")
        return pd.DataFrame(rows)