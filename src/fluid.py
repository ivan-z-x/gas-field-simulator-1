from __future__ import annotations
import math
import pandas as pd
from src.interpolator import LinearInterpolator


class Fluid:
    """PVT-модель природного газа. Z-фактор по GERG-91 (ГОСТ 30319.2-96)."""

    R     = 8.314462618
    P_STD = 101325.0
    T_STD = 293.15
    _A = [0.3265, -1.0700, -0.5339, 0.01569, -0.05165,
          0.5475, -0.7361,  0.1844,  0.1056,  0.6134,  0.7210]

    def __init__(self, M: float, rho_c: float, xa: float, xy: float,
                 T: float, interp_csv_path: str):
        self.M = M; self.rho_c = rho_c; self.xa = xa; self.xy = xy; self.T = T

        g = float(rho_c); y_n2 = float(xa); y_co2 = float(xy)
        y_hc = max(0.0, 1.0 - y_n2 - y_co2)
        Tpc_hc = (168.0 + 325.0*g - 12.5*g**2) * 5.0/9.0
        Ppc_hc = (677.0 + 15.0*g  - 37.5*g**2) / 14.6959
        Tc = {"N2": 126.2, "CO2": 304.2}; Pc = {"N2": 33.98, "CO2": 72.8}
        self.Tpc = y_n2*Tc["N2"] + y_co2*Tc["CO2"] + y_hc*Tpc_hc
        self.Ppc = y_n2*Pc["N2"] + y_co2*Pc["CO2"] + y_hc*Ppc_hc

        # Читаем CSV: автоматически определяем разделитель (; или ,)
        try:
            df = pd.read_csv(interp_csv_path, sep=';')
            if df.shape[1] < 2:
                df = pd.read_csv(interp_csv_path, sep=',')
        except Exception:
            df = pd.read_csv(interp_csv_path, sep=',')

        cols = df.columns.tolist()
        self._mu_interp = LinearInterpolator(
            xs=df[cols[0]].astype(float).tolist(),
            ys=df[cols[1]].astype(float).tolist(),
        )

    def _z_from_rho(self, rho: float, Tpr: float) -> float:
        A = self._A
        c1 = A[0]+A[1]/Tpr+A[2]/Tpr**3+A[3]/Tpr**4+A[4]/Tpr**5
        c2 = A[5]+A[6]/Tpr+A[7]/Tpr**2
        c3 = A[8]*(A[5]+A[6]/Tpr+A[7]/Tpr**2)
        e  = math.exp(-A[10]*rho**2)
        c4 = A[9]*(1.0+A[10]*rho**2)*(rho**2/Tpr**3)*e
        return 1.0 + c1*rho + c2*rho**2 - c3*rho**5 + c4

    def z(self, P: float) -> float:
        P = max(float(P), 1e-6)
        Ppr = P/self.Ppc; Tpr = self.T/self.Tpc
        def eq(rho):
            if rho < 1e-12: return 1.0
            return self._z_from_rho(rho, Tpr) - 0.27*Ppr/(Tpr*rho)
        lo, hi = 1e-6, 0.99
        if eq(lo)*eq(hi) > 0:
            for t in [0.1,0.2,0.3,0.5,0.7,0.9]:
                if eq(lo)*eq(t) < 0: hi=t; break
            else: return 0.85
        for _ in range(80):
            mid = (lo+hi)/2.0
            if eq(lo)*eq(mid) < 0: hi=mid
            else: lo=mid
            if hi-lo < 1e-10: break
        return max(0.2, min(2.0, self._z_from_rho((lo+hi)/2.0, Tpr)))

    def ro(self, P: float) -> float:
        return float(P)*101325.0*self.M/(self.z(P)*self.R*self.T)

    def ro_std(self) -> float:
        return self.P_STD*self.M/(self.R*self.T_STD)

    def bg(self, P: float) -> float:
        return (self.P_STD*self.z(P)*self.T)/(float(P)*101325.0*self.T_STD)

    def mu(self, P: float) -> float:
        return self._mu_interp.predict(float(P))
