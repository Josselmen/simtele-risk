"""SIMTELE-RISK: integracion (trapecio, Simpson, Romberg). Proyecto educativo, datos sinteticos."""
import math
import numpy as np

NAN = float("nan")

# --------------------------------------------------------------- INTEGRACION
def trap_arr(y, h):
    y = np.asarray(y, float)
    return float(h * (0.5 * y[0] + y[1:-1].sum() + 0.5 * y[-1]))


def simpson_arr(y, h):
    y = np.asarray(y, float); n = len(y) - 1
    if n % 2:
        raise ValueError("Simpson 1/3 compuesto requiere n par")
    return float(h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum()))


def romberg_arr(y, h):
    """y con 2^m+1 puntos equiespaciados. Devuelve la tabla de Romberg R[i, j]."""
    y = np.asarray(y, float); N = len(y) - 1; m = int(round(math.log2(N)))
    if 2 ** m != N:
        raise ValueError("Romberg requiere 2^m subintervalos")
    R = np.zeros((m + 1, m + 1))
    for i in range(m + 1):
        s = 2 ** (m - i)
        R[i, 0] = trap_arr(y[::s], h * s)
        for j in range(1, i + 1):
            R[i, j] = R[i, j - 1] + (R[i, j - 1] - R[i - 1, j - 1]) / (4 ** j - 1)
    return R


def grid(f, a, b, n):
    x = np.linspace(a, b, n + 1)
    return x, np.array([f(t) for t in x]), (b - a) / n


