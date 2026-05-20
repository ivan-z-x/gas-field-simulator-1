# Линейная интерполяция без NumPy (только встроенные структуры Python)
class LinearInterpolator:
    def __init__(self, xs: list[float], ys: list[float]):
        if len(xs) != len(ys):
            raise ValueError("xs и ys должны иметь одинаковую длину")
        if len(xs) < 2:
            raise ValueError("Нужно минимум две точки")
        self.xs = list(xs)
        self.ys = list(ys)
        # Проверка строгого возрастания xs
        for i in range(len(self.xs) - 1):
            if self.xs[i] >= self.xs[i + 1]:
                raise ValueError("xs должны строго возрастать")

    def predict(self, xp: float) -> float:
        # Экстраполяция запрещена (только интерполяция)
        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError(f"xp={xp} вне диапазона [{self.xs[0]}, {self.xs[-1]}]")
        if xp == self.xs[-1]:
            return self.ys[-1]
        for i in range(len(self.xs) - 1):
            x0, x1 = self.xs[i], self.xs[i + 1]
            if x0 <= xp <= x1:
                y0, y1 = self.ys[i], self.ys[i + 1]
                return y0 + (y1 - y0) * (xp - x0) / (x1 - x0)
        raise RuntimeError("Не найден интервал для интерполяции")