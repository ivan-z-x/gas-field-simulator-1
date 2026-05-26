class DCS:
    """Дожимная компрессорная станция (ДКС)."""

    def __init__(self, CR: float, P_line: float, q_ext: float = 0.0):
        """
        CR     — степень сжатия (≥ 1.0)
        P_line — давление в магистрали [атм]
        q_ext  — сторонний газ на манифолде [ст.м³/сут]
        """
        if CR < 1.0:
            raise ValueError("CR должен быть >= 1.0")
        self.CR = CR
        self.P_line = P_line
        self.q_ext = q_ext

    def P_in(self) -> float:
        """Давление на входе в ДКС [атм]."""
        return self.P_line / self.CR
