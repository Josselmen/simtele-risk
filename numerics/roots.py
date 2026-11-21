"""SIMTELE-RISK: raices (biseccion, Newton, secante). Proyecto educativo, datos sinteticos."""
import math
import numpy as np

NAN = float("nan")

# -------------------------------------------------------------------- RAICES
def bisection(f, a, b, tol=1e-10, maxit=200, xref=None):
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("f(a) y f(b) tienen el mismo signo: no hay garantia de raiz")
    hist = []
    for k in range(1, maxit + 1):
        c = 0.5 * (a + b); fc = f(c)
        hist.append(dict(k=k, a=a, b=b, x=c, fx=fc, step=0.5 * (b - a),
                         err=(abs(c - xref) if xref is not None else NAN)))
        if abs(fc) < tol or 0.5 * (b - a) < tol:
            return dict(root=c, iters=hist, converged=True)
        if fa * fc < 0: b, fb = c, fc
        else: a, fa = c, fc
    return dict(root=c, iters=hist, converged=False)


def newton(f, df, x0, tol=1e-10, maxit=100, xref=None, m=1.0, use_res=True):
    x = x0; hist = []
    for k in range(1, maxit + 1):
        d = df(x)
        if abs(d) < 1e-14:
            return dict(root=x, iters=hist, converged=False, reason="derivada ~ 0")
        xn = x - m * f(x) / d
        hist.append(dict(k=k, x=xn, fx=f(xn), step=abs(xn - x),
                         err=(abs(xn - xref) if xref is not None else NAN)))
        done = abs(xn - x) < tol or (use_res and abs(f(xn)) < tol)
        x = xn
        if done:
            return dict(root=x, iters=hist, converged=True, reason="ok")
    return dict(root=x, iters=hist, converged=False, reason="maxit")


def secant(f, x0, x1, tol=1e-10, maxit=100, xref=None, use_res=True):
    hist = []
    for k in range(1, maxit + 1):
        f0, f1 = f(x0), f(x1)
        if abs(f1 - f0) < 1e-300:
            return dict(root=x1, iters=hist, converged=False, reason="f1 = f0")
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        hist.append(dict(k=k, x=x2, fx=f(x2), step=abs(x2 - x1),
                         err=(abs(x2 - xref) if xref is not None else NAN)))
        done = abs(x2 - x1) < tol or (use_res and abs(f(x2)) < tol)
        x0, x1 = x1, x2
        if done:
            return dict(root=x1, iters=hist, converged=True, reason="ok")
    return dict(root=x1, iters=hist, converged=False, reason="maxit")


