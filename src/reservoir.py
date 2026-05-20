from __future__ import annotations
from dataclasses import dataclass
from src.fluid import Fluid

@dataclass
class ResProps:
    P: float      # атм (начальное давление)
    V: float      # м³ (объём пласта – поровый объём, но здесь используется для GIIP)
    T: float      # K

class Reservoir:
    def __init__(self, resprops: ResProps, fluid: Fluid):
        self.resprops = resprops
        self.fluid = fluid
        # Вычисление начальных геологических запасов газа (GIIP) в стандартных м³
        Pi = resprops.P
        zi = self.fluid.z(Pi)
        T_std = fluid.T_STD_K
        P_std = fluid.P_STD_PA / 101325.0   # атм
        self.GIIP = resprops.V * (Pi / P_std) * (T_std / resprops.T) * (1.0 / zi)

        # Накопленная добыча (инициализация)
        self.Gp = 0.0

    def p2(self, q_total: float, dt: float = 1.0) -> float:
        """
        Обновление пластового давления по материальному балансу:
        P/z = (P_i / z_i) * (1 - Gp / GIIP)
        q_total – дебит [ст.м³/сут], dt – шаг [сут].
        """
        # Обновляем накопленную добычу
        self.Gp += q_total * dt

        if self.Gp >= self.GIIP:
            return 1.0

        P_i = self.resprops.P
        z_i = self.fluid.z(P_i)
        P_z_initial = P_i / z_i
        P_z_new = P_z_initial * (1.0 - self.Gp / self.GIIP)

        # Итерационное решение P = P_z_new * z(P)
        P_new = P_i
        for _ in range(10):
            z_new = self.fluid.z(P_new)
            f = P_new - P_z_new * z_new
            if abs(f) < 1e-4:
                break
            P_new = P_z_new * z_new   # простая итерация
        return max(1.0, P_new)