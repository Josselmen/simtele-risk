"""Semana 1: pruebas de teoria de errores y sistemas lineales contra soluciones conocidas."""
import math

import numpy as np
import pytest

import numerics as nm


def test_errores():
    assert nm.err_abs(1.41, math.sqrt(2)) == pytest.approx(0.004213562, rel=1e-5)
    assert math.isnan(nm.err_rel(1.0, 0.0))          # no se divide por cero
    assert nm.err_pct(110, 100) == pytest.approx(10.0)


def test_residuo_y_amplificacion():
    A = np.array([[2.0, 0], [0, 4.0]]); b = np.array([2.0, 4.0])
    assert nm.residual(A, b, np.array([1.0, 1.0])) == 0.0
    assert nm.amplification(10, 11, 100, 110) == pytest.approx(1.0)


def test_sistemas_lineales_mismo_resultado():
    A = nm.tridiag(4, 4.0, 1.0); xs = np.array([1, 2, 3, 4.0]); b = A @ xs
    np.testing.assert_allclose(b, [6, 12, 18, 19])
    for r in (nm.jacobi(A, b, xref=xs), nm.gauss_seidel(A, b, xref=xs), nm.cg(A, b, xref=xs)):
        assert r["converged"]
        np.testing.assert_allclose(r["x"], xs, atol=1e-8)
    np.testing.assert_allclose(nm.lu_solve(A, b), xs, atol=1e-12)


def test_lu_reproduce_PA():
    A = np.array([[1.0, 2, 3], [4, 5, 6], [7, 8, 10]])
    L, U, perm = nm.lu_decomp(A)
    np.testing.assert_allclose(L @ U, A[perm], atol=1e-12)


def test_cg_termina_en_n_pasos_y_gs_mas_rapido_que_jacobi():
    A = nm.tridiag(4, 4.0, 1.0); b = A @ np.arange(1, 5.0)
    assert nm.cg(A, b)["n"] == 4
    assert nm.gauss_seidel(A, b)["n"] < nm.jacobi(A, b)["n"]


def test_hilbert_jacobi_diverge():
    H = nm.hilbert(4); b = H @ np.ones(4)
    assert not nm.jacobi(H, b, maxit=2000)["converged"]
    assert np.linalg.cond(H) == pytest.approx(1.55e4, rel=0.01)


def test_cg_rechaza_matriz_no_spd():
    B = np.array([[1, 3, 0, 0], [0, 1, 2, 0], [5, 0, 1, 1], [0, 2, 0, 1.0]])
    assert not nm.cg(B, B @ np.arange(1, 5.0), maxit=40)["converged"]
