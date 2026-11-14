"""SIMTELE-RISK: ecuaciones diferenciales (RK2, RK4). Proyecto educativo, datos sinteticos."""
import math
import numpy as np

NAN = float("nan")

# ----------------------------------------------------------------------- EDO
def _march(step, f, t0, y0, h, tend):
    t = t0; y = np.array(y0, float); T = [t]; Y = [y.copy()]
    while t < tend - 1e-12:
        hh = min(h, tend - t)
        y = step(f, t, y, hh); t += hh
        T.append(t); Y.append(y.copy())
    return np.array(T), np.array(Y)


def _rk2(f, t, y, h):  # RK2 variante PUNTO MEDIO
    k1 = f(t, y); k2 = f(t + h / 2, y + h / 2 * k1)
    return y + h * k2


def _rk4(f, t, y, h):
    k1 = f(t, y); k2 = f(t + h / 2, y + h / 2 * k1)
    k3 = f(t + h / 2, y + h / 2 * k2); k4 = f(t + h, y + h * k3)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def rk2(f, t0, y0, h, tend): return _march(_rk2, f, t0, y0, h, tend)
def rk4(f, t0, y0, h, tend): return _march(_rk4, f, t0, y0, h, tend)
