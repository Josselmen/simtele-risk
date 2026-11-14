"""SIMTELE-RISK: teoria de errores. Proyecto educativo, datos sinteticos."""
import math
import numpy as np

NAN = float("nan")

# ------------------------------------------------------------------ ERRORES
def err_abs(x, xr):
    return float(np.linalg.norm(np.atleast_1d(x) - np.atleast_1d(xr), np.inf))


def err_rel(x, xr):
    """Error relativo; devuelve NaN si |xr| ~ 0 (evita dividir por cero)."""
    n = float(np.linalg.norm(np.atleast_1d(xr), np.inf))
    return NAN if n < 1e-14 else err_abs(x, xr) / n


def err_pct(x, xr):
    return 100.0 * err_rel(x, xr)


def residual(A, b, x):
    return float(np.linalg.norm(A @ x - b, np.inf))


def amplification(out_ref, out_pert, in_ref, in_pert):
    """Factor de amplificacion relativo = (dSalida/Salida) / (dEntrada/Entrada)."""
    return (abs(out_pert - out_ref) / abs(out_ref)) / (abs(in_pert - in_ref) / abs(in_ref))


