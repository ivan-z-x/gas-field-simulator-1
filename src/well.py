# src/well.py
# Модель скважины: приток (IPR)

import math
from src.fluid import Fluid

class Well:
    def __init__(self, fluid: Fluid, k: float, h: float, re: float, rw: float):
        self.fluid = fluid
        self.k = k          # Проницаемость, мД
        self.h = h          # Эффективная мощность, м
        self.re = re        # Радиус контура питания, м
        self.rw = rw        # Радиус скважины, м
        # Коэффициент перевода единиц β = 0.00852702 ст.м³/(сут·мД·м·атм)
        self.beta = 0.00852702

    def ipr(self, P_res: float, P_bhp: float) -> float:
        """
        Дебит скважины по закону Дарси (радиальный приток).
        P_res - пластовое давление [атм]
        P_bhp - забойное давление [атм]
        Возвращает дебит [ст.м³/сут] (положительный, если P_res > P_bhp, иначе 0)
        """
        if P_bhp >= P_res:
            return 0.0

        # Вязкость при пластовом давлении
        mu = self.fluid.mu(P_res)  # сП
        # Коэффициент продуктивности
        C = self.beta * self.k * self.h / (mu * math.log(self.re / self.rw))
        # Формула притока
        q = C * (P_res - P_bhp)
        return q