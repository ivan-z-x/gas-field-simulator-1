from dataclasses import dataclass

@dataclass
class NodeState:
    name: str
    P_in: float
    P_out: float
    dP: float
    q_std: float
    q_res: float | None = None
    v: float | None = None
    rho: float | None = None