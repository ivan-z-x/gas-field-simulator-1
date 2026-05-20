from __future__ import annotations

import math
import pandas as pd
from src.interpolator import LinearInterpolator

class Fluid:
    """
    PVT-модель природного газа.
    Реализует модифицированное уравнение GERG-91 для z-фактора.
    Индивидуальные данные студента Али Иван: xa=0.3745, xy=0.9507, rho_c=0.6799.
    """
    R = 8.314462618  # Дж/(моль·К)
    P_STD_PA = 101325.0   # Па
    T_STD_K = 293.15      # К

    def __init__(self, interp_csv_path: str):
        # Индивидуальные параметры (Али Иван) – жёстко заданы
        self.M = 0.0150                     # кг/моль (молярная масса природного газа)
        self.rho_c = 0.6799                 # относительная плотность по воздуху
        self.xa = 0.3745                    # мольная доля N2, %
        self.xy = 0.9507                    # мольная доля CO2, %
        self.T = 310.0                      # К (пластовая температура)

        # Чтение данных для интерполяции вязкости
        df = pd.read_csv(interp_csv_path)
        self._mu_interp = LinearInterpolator(
            xs=df["pressure, atm"].astype(float).tolist(),
            ys=df["viscosity, cP"].astype(float).tolist(),
        )

        # Псевдокритические свойства (правило смешения Кея)
        y_n2 = self.xa / 100.0
        y_co2 = self.xy / 100.0
        y_hc = max(0.0, 1.0 - y_n2 - y_co2)

        # Псевдокритические параметры углеводородов (корреляция Саттона)
        gamma_g = self.rho_c
        Tpc_hc_R = 168.0 + 325.0 * gamma_g - 12.5 * gamma_g**2   # Ранкин
        Ppc_hc_psia = 677.0 + 15.0 * gamma_g - 37.5 * gamma_g**2  # psia
        Tpc_hc_K = Tpc_hc_R * 5.0 / 9.0
        Ppc_hc_atm = Ppc_hc_psia / 14.6959

        # Критические параметры N2 и CO2
        T_crit = {"N2": 126.2, "CO2": 304.2, "HC": Tpc_hc_K}
        P_crit = {"N2": 33.98, "CO2": 72.8, "HC": Ppc_hc_atm}

        self.Tpc = y_n2 * T_crit["N2"] + y_co2 * T_crit["CO2"] + y_hc * T_crit["HC"]
        self.Ppc = y_n2 * P_crit["N2"] + y_co2 * P_crit["CO2"] + y_hc * P_crit["HC"]

    def z(self, P: float) -> float:
        """
        Коэффициент сжимаемости по модифицированному GERG-91 (корреляция Дранчука-Абу-Кассема).
        P в атмосферах.
        """
        P = max(float(P), 1e-6)
        Ppr = P / self.Ppc
        Tpr = self.T / self.Tpc

        # Коэффициенты DAK
        A1 = 0.3265; A2 = -1.0700; A3 = -0.5339; A4 = 0.01569; A5 = -0.05165
        A6 = 0.5475; A7 = -0.7361; A8 = 0.1844; A9 = 0.1056; A10 = 0.6134; A11 = 0.7210

        rho_r = 0.25  # начальное приближение приведённой плотности
        for _ in range(20):
            term1 = A1 + A2/Tpr + A3/Tpr**3 + A4/Tpr**4 + A5/Tpr**5
            term2 = (A6 + A7/Tpr + A8/Tpr**2) * rho_r
            term3 = (A9 + A10/Tpr) * rho_r**2
            term4 = A11 * rho_r**5 / Tpr**3
            f = rho_r + term2 + term3 + term4 - Ppr / (term1 * rho_r)
            f_prime = 1 + 2*term2 + 3*term3 + 6*term4 - (Ppr/term1)*(-1/rho_r**2)
            rho_r_new = rho_r - f / f_prime
            if abs(rho_r_new - rho_r) < 1e-8:
                rho_r = rho_r_new
                break
            rho_r = rho_r_new
        z = Ppr / (rho_r * term1)
        return max(0.5, min(1.5, z))

    def ro(self, P: float) -> float:
        """Плотность реального газа [кг/м³], P в атмосферах."""
        P_pa = P * 101325.0
        return P_pa * self.M / (self.z(P) * self.R * self.T)

    def ro_ideal(self, P: float) -> float:
        """Плотность идеального газа [кг/м³] для сравнения."""
        P_pa = P * 101325.0
        return P_pa * self.M / (self.R * self.T)

    def bg(self, P: float) -> float:
        """Объёмный коэффициент газа [ст.м³/пл.м³]."""
        return (self.P_STD_PA * self.z(P) * self.T) / (P * 101325.0 * self.T_STD_K)

    def mu(self, P: float) -> float:
        """Вязкость [сП] – линейная интерполяция табличных данных."""
        return self._mu_interp.predict(float(P))