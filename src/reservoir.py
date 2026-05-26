from __future__ import annotations
from dataclasses import dataclass
from src.fluid import Fluid


@dataclass
class ResProps:
    P: float    # текущее пластовое давление [атм]
    V: float    # поровый объём пласта [м³]
    T: float    # пластовая температура [К]


class Reservoir:
    """Модель пласта: материальный баланс."""

    def __init__(self, resprops: ResProps, fluid: Fluid):
        self.resprops = resprops
        self.fluid = fluid

    def p2(self, q_total: float, dt: float = 1.0) -> float:
        """
        Новое пластовое давление после шага dt [сут].
        q_total — суммарный дебит скважин куста [ст.м³/сут] (без q_ext).
        Формула: P_new = P - Z*rho_std/rho * q_total/V * dt
        resprops.P не изменяется здесь — изменяет FieldSimulator.
        """
        P = self.resprops.P
        z_p   = self.fluid.z(P)
        rho_p = self.fluid.ro(P)             # кг/м³ при пластовых условиях
        rho_s = self.fluid.ro_std()           # кг/м³ при стандартных условиях

        # q_total [ст.м³/сут] → [м³/с] не нужен: dt в сутках, V в м³
        dP = z_p * rho_s / rho_p * q_total / self.resprops.V * dt
        P_new = P - dP
        return max(0.1, P_new)
