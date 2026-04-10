class LinearInterpolator:
    def __init__(self, xs: list, ys: list):
        if len(xs) != len(ys):
            raise ValueError("xs и ys должны быть одинаковой длины")
        if not all(xs[i] < xs[i+1] for i in range(len(xs)-1)):
            raise ValueError("xs должен быть строго возрастающим")
        self.xs = xs
        self.ys = ys

    def predict(self, xp: float) -> float:
        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError(f"xp={xp} вне диапазона [{self.xs[0]}, {self.xs[-1]}]")
        lo, hi = 0, len(self.xs)-1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if xp < self.xs[mid]:
                hi = mid
            else:
                lo = mid
        x1, y1 = self.xs[lo], self.ys[lo]
        x2, y2 = self.xs[hi], self.ys[hi]
        return y1 + (y2 - y1) * (xp - x1) / (x2 - x1)