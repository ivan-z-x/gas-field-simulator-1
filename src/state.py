from dataclasses import dataclass


@dataclass
class NodeState:
    name: str
    P_in: float        # давление на входе [атм]
    P_out: float       # давление на выходе [атм]
    dP: float          # перепад давления [атм]
    q_std: float       # дебит при стандартных условиях [ст.м³/сут]
    q_res: float | None = None   # дебит при пластовых условиях [м³/сут]
    v: float | None = None       # скорость потока [м/с]
    rho: float | None = None     # плотность газа [кг/м³]
