# src/fluid.py
# Свойства реального газа

import csv
from src.interpolator import LinearInterpolator

class Fluid:
    def __init__(self, M: float, rho_c: float, xa: float, xy: float, T: float):
        # Параметры газа (Али Иван)
        self.M = M
        self.rho_c = rho_c
        self.xa = xa
        self.xy = xy
        self.T = T

        # Загрузка таблицы вязкости из interp_data.csv
        pressures, viscosities = [], []
        with open('interp_data.csv', 'r') as f:
            reader = csv.reader(f, delimiter=';')
            next(reader)  # пропускаем заголовок
            for row in reader:
                p = float(row[0].replace(',', '.'))
                mu = float(row[1].replace(',', '.'))
                pressures.append(p)
                viscosities.append(mu)
        self.mu_interp = LinearInterpolator(pressures, viscosities)

    def mu(self, P: float) -> float:
        """Динамическая вязкость [сП] при давлении P (атм)"""
        return self.mu_interp.predict(P)

    # Временные заглушки 
    def z(self, P: float) -> float:
        return 1.0

    def ro(self, P: float) -> float:
        return 1.0

    def bg(self, P: float) -> float:
        return 1.0