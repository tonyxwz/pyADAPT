import numpy as np
from multiprocessing import Pool


class TestMpEval(object):
    def __init__(self):
        self.text_formula = "sin(x) * cos(x)"
        self.formula = compile(self.text_formula, "<string>", "eval")

    def compute(self, x):
        return eval(self.formula, {"sin": np.sin, "cos": np.cos, "x": x})


if __name__ == "__main__":
    pool = Pool(4)
    c = TestMpEval()
    res = list()
    t = np.linspace(-2 * np.pi, 2 * np.pi, 100)
    for i in t:
        res.append(pool.apply_async(c.compute, args=(i)))
    y = np.zeros(100)
    for i, r in enumerate(res):
        y[i] = r.get()
    import matplotlib.pyplot as plt
    plt.plot(t, y)
    plt.show()