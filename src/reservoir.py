# src/reservoir.py
# Модель пласта (материальный баланс)

from dataclasses import dataclass
from src.fluid import Fluid

@dataclass
class ResProps:
    """Контейнер свойств пласта."""
    P: float   # Давление, атм
    V: float   # Объем порового пространства, м³
    T: float   # Температура, К

class Reservoir:
    def __init__(self, resprops: ResProps, fluid: Fluid):
        self.resprops = resprops
        self.fluid = fluid
        # Стандартная плотность (при Pstd, Tstd)
        self.rho_std = fluid.ro(1.0)

    def p2(self, q_total_std: float, dt: float = 1.0) -> float:
        """
        Рассчитывает новое пластовое давление через dt суток.
        q_total_std - суммарный дебит скважин [ст.м³/сут]
        Возвращает P_res на следующий шаг.
        """
        P = self.resprops.P
        V_res = self.resprops.V
        # Плотность и Z-фактор при текущем P
        rho_res = self.fluid.ro(P)
        Z = self.fluid.z(P)

        # Масса газа, добытая за dt
        mass_produced = q_total_std * dt * self.rho_std
        # Начальная масса в пласте
        mass_res = rho_res * V_res
        # Новая масса
        mass_new = max(mass_res - mass_produced, 0.0)

        # Новое давление из уравнения состояния: P_new = (mass_new * Z_new * R * T) / (M * V_res)
        # Начальное приближение, используя старый Z
        P_new = mass_new * Z * self.fluid.R * self.resprops.T / (self.fluid.M * V_res) / 101325.0

        # Уточнение: пара итераций для учета изменения Z от давления
        for _ in range(3):
            Z_new = self.fluid.z(P_new)
            P_new = mass_new * Z_new * self.fluid.R * self.resprops.T / (self.fluid.M * V_res) / 101325.0
        return max(P_new, 0.0)