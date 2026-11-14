"""SIMTELE-RISK: sistemas lineales (LU, Jacobi, Gauss-Seidel, gradiente conjugado). Proyecto educativo, datos sinteticos."""
import math
import numpy as np

NAN = float("nan")

from .errors import err_abs, residual

# --------------------------------------------------------- SISTEMAS LINEALES
def lu_decomp(A):
    """LU con pivoteo parcial: P*A = L*U. Devuelve L, U, perm."""
    A = np.array(A, dtype=float)
    n = len(A)
    perm = np.arange(n)
    for k in range(n - 1):
        p = k + int(np.argmax(np.abs(A[k:, k])))
        if abs(A[p, k]) < 1e-300:
            raise ValueError("Matriz singular")
        if p != k:
            A[[k, p]] = A[[p, k]]
            perm[[k, p]] = perm[[p, k]]
        for i in range(k + 1, n):
            A[i, k] /= A[k, k]
            A[i, k + 1:] -= A[i, k] * A[k, k + 1:]
    return np.tril(A, -1) + np.eye(n), np.triu(A), perm


def lu_solve(A, b):
    L, U, perm = lu_decomp(A)
    b = np.array(b, dtype=float)[perm]
    n = len(b)
    y = np.zeros(n)
    for i in range(n):
        y[i] = b[i] - L[i, :i] @ y[:i]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - U[i, i + 1:] @ x[i + 1:]) / U[i, i]
    return x


def _rec(k, x, A, b, xref):
    return dict(k=k, x=x.copy(), res=residual(A, b, x),
                err=(err_abs(x, xref) if xref is not None else NAN))


def jacobi(A, b, tol=1e-10, maxit=5000, xref=None):
    A = np.array(A, float); b = np.array(b, float)
    n = len(b); x = np.zeros(n); D = np.diag(A); R = A - np.diagflat(D)
    nb = np.linalg.norm(b, np.inf); hist = []
    for k in range(1, maxit + 1):
        x = (b - R @ x) / D
        r = _rec(k, x, A, b, xref); hist.append(r)
        if not np.isfinite(r["res"]) or r["res"] > 1e100:
            return dict(x=x, iters=hist, converged=False, n=k, reason="diverge")
        if r["res"] / nb < tol:
            return dict(x=x, iters=hist, converged=True, n=k, reason="ok")
    return dict(x=x, iters=hist, converged=False, n=maxit, reason="maxit")


def gauss_seidel(A, b, tol=1e-10, maxit=5000, xref=None):
    A = np.array(A, float); b = np.array(b, float)
    n = len(b); x = np.zeros(n); nb = np.linalg.norm(b, np.inf); hist = []
    for k in range(1, maxit + 1):
        for i in range(n):
            x[i] = (b[i] - A[i, :i] @ x[:i] - A[i, i + 1:] @ x[i + 1:]) / A[i, i]
        r = _rec(k, x, A, b, xref); hist.append(r)
        if not np.isfinite(r["res"]) or r["res"] > 1e100:
            return dict(x=x, iters=hist, converged=False, n=k, reason="diverge")
        if r["res"] / nb < tol:
            return dict(x=x, iters=hist, converged=True, n=k, reason="ok")
    return dict(x=x, iters=hist, converged=False, n=maxit, reason="maxit")


def cg(A, b, tol=1e-10, maxit=None, xref=None):
    """Gradiente conjugado estandar (exige A simetrica definida positiva)."""
    A = np.array(A, float); b = np.array(b, float)
    n = len(b); maxit = maxit or 10 * n
    x = np.zeros(n); r = b.copy(); p = r.copy(); rs = r @ r
    nb = np.linalg.norm(b, np.inf); hist = []
    for k in range(1, maxit + 1):
        Ap = A @ p; pAp = p @ Ap
        if pAp <= 0 or not np.isfinite(pAp):
            return dict(x=x, iters=hist, converged=False, n=k - 1, reason="p'Ap<=0 (A no SPD)")
        alpha = rs / pAp
        x = x + alpha * p; r = r - alpha * Ap; rsn = r @ r
        rec = _rec(k, x, A, b, xref); hist.append(rec)
        if rec["res"] / nb < tol:
            return dict(x=x, iters=hist, converged=True, n=k, reason="ok")
        p = r + (rsn / rs) * p; rs = rsn
    return dict(x=x, iters=hist, converged=False, n=maxit, reason="maxit")


def hilbert(n):
    return np.array([[1.0 / (i + j + 1) for j in range(n)] for i in range(n)])


def tridiag(n, d, o):
    return d * np.eye(n) + o * (np.eye(n, k=1) + np.eye(n, k=-1))


