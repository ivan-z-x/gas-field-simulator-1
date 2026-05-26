class LinearInterpolator:
    """Линейная интерполяция (с экстраполяцией крайними значениями)."""

    def __init__(self, xs: list[float], ys: list[float]):
        if len(xs) != len(ys) or len(xs) < 2:
            raise ValueError("Нужно минимум 2 точки с одинаковой длиной xs и ys")
        for i in range(len(xs) - 1):
            if xs[i] >= xs[i + 1]:
                raise ValueError("xs должны строго возрастать")
        self.xs = list(xs)
        self.ys = list(ys)

    def predict(self, xp: float) -> float:
        # Экстраполяция крайними значениями вместо исключения
        if xp <= self.xs[0]:
            return self.ys[0]
        if xp >= self.xs[-1]:
            return self.ys[-1]
        for i in range(len(self.xs) - 1):
            x0, x1 = self.xs[i], self.xs[i + 1]
            if x0 <= xp <= x1:
                y0, y1 = self.ys[i], self.ys[i + 1]
                return y0 + (y1 - y0) * (xp - x0) / (x1 - x0)
        return self.ys[-1]
