class DCS:
    def __init__(self, CR: float, P_line: float, q_ext: float = 0.0):
        if CR < 1.0:
            raise ValueError("CR должен быть >= 1.0")
        self.CR = CR          # степень сжатия
        self.P_line = P_line  # давление на выходе (линия)
        self.q_ext = q_ext    # внешний расход (например, от соседнего куста)

    def P_in(self) -> float:
        """Давление на входе в ДКС."""
        return self.P_line / self.CR