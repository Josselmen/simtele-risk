"""SIMTELE-RISK: interpolacion (Lagrange, Newton, spline cubico natural). Proyecto educativo, datos sinteticos."""
import math
import numpy as np

NAN = float("nan")

from .linear import lu_solve

# ------------------------------------------------------------- INTERPOLACION
def lagrange(xn, yn, x):
    xn = np.asarray(xn, float); yn = np.asarray(yn, float)
    x = np.atleast_1d(np.asarray(x, float)); out = np.zeros_like(x)
    for i in range(len(xn)):
        Li = np.ones_like(x)
        for j in range(len(xn)):
            if j != i:
                Li *= (x - xn[j]) / (xn[i] - xn[j])
        out += yn[i] * Li
    return out


def newton_table(xn, yn):
    """Tabla de diferencias divididas; los coeficientes son T[0, :]."""
    xn = np.asarray(xn, float); n = len(xn)
    T = np.zeros((n, n)); T[:, 0] = yn
    for j in range(1, n):
        for i in range(n - j):
            T[i, j] = (T[i + 1, j - 1] - T[i, j - 1]) / (xn[i + j] - xn[i])
    return T


def newton_eval(xn, T, x):
    xn = np.asarray(xn, float); c = T[0, :]
    x = np.atleast_1d(np.asarray(x, float)); n = len(xn)
    p = np.full_like(x, c[-1])
    for k in range(n - 2, -1, -1):
        p = p * (x - xn[k]) + c[k]
    return p


def spline_natural(xn, yn):
    """Spline cubico natural (S''=0 en los extremos). Devuelve (evaluador, M)."""
    xn = np.asarray(xn, float); yn = np.asarray(yn, float); n = len(xn)
    h = np.diff(xn); M = np.zeros(n)
    if n > 2:
        A = np.zeros((n - 2, n - 2)); rhs = np.zeros(n - 2)
        for i in range(1, n - 1):
            A[i - 1, i - 1] = 2 * (h[i - 1] + h[i])
            if i > 1: A[i - 1, i - 2] = h[i - 1]
            if i < n - 2: A[i - 1, i] = h[i]
            rhs[i - 1] = 6 * ((yn[i + 1] - yn[i]) / h[i] - (yn[i] - yn[i - 1]) / h[i - 1])
        M[1:-1] = lu_solve(A, rhs)

    def S(x):
        x = np.atleast_1d(np.asarray(x, float))
        i = np.clip(np.searchsorted(xn, x) - 1, 0, n - 2)
        hi = h[i]; a = xn[i + 1] - x; b = x - xn[i]
        return (M[i] * a ** 3 / (6 * hi) + M[i + 1] * b ** 3 / (6 * hi)
                + (yn[i] / hi - M[i] * hi / 6) * a + (yn[i + 1] / hi - M[i + 1] * hi / 6) * b)
    return S, M


def lebesgue(xn, grid=4001):
    xn = np.asarray(xn, float); xs = np.linspace(xn[0], xn[-1], grid); s = np.zeros_like(xs)
    for i in range(len(xn)):
        Li = np.ones_like(xs)
        for j in range(len(xn)):
            if j != i:
                Li *= (xs - xn[j]) / (xn[i] - xn[j])
        s += np.abs(Li)
    return float(s.max())


