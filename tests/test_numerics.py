"""Pruebas de los metodos numericos contra soluciones analiticas conocidas."""
import math

import numpy as np
import pytest

import numerics as nm


def test_errores():
    assert nm.err_abs(1.41, math.sqrt(2)) == pytest.approx(0.004213562, rel=1e-5)
    assert math.isnan(nm.err_rel(1.0, 0.0))          # no se divide por cero
    assert nm.err_pct(110, 100) == pytest.approx(10.0)


def test_sistemas_lineales_mismo_resultado():
    A = nm.tridiag(4, 4.0, 1.0); xs = np.array([1, 2, 3, 4.0]); b = A @ xs
    np.testing.assert_allclose(b, [6, 12, 18, 19])
    for r in (nm.jacobi(A, b, xref=xs), nm.gauss_seidel(A, b, xref=xs), nm.cg(A, b, xref=xs)):
        assert r["converged"]
        np.testing.assert_allclose(r["x"], xs, atol=1e-8)
    np.testing.assert_allclose(nm.lu_solve(A, b), xs, atol=1e-12)


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


def test_raices_x2_menos_2():
    f = lambda x: x * x - 2; df = lambda x: 2 * x
    for r in (nm.bisection(f, 1, 2, 1e-12), nm.newton(f, df, 1.5, 1e-12), nm.secant(f, 1, 2, 1e-12)):
        assert r["converged"] and r["root"] == pytest.approx(math.sqrt(2), abs=1e-10)
    with pytest.raises(ValueError):
        nm.bisection(f, 2, 3)                          # sin cambio de signo


def test_integracion_seno():
    for n, tol in ((64, 1e-2), (64, 1e-5)):
        x, y, h = nm.grid(np.sin, 0, math.pi, n)
    assert abs(nm.trap_arr(y, h) - 2) < 1e-2
    assert abs(nm.simpson_arr(y, h) - 2) < 1e-6
    assert abs(nm.romberg_arr(y, h)[-1, -1] - 2) < 1e-12
    with pytest.raises(ValueError):
        nm.simpson_arr(np.ones(4), 1.0)                # n impar


def test_ordenes_edo():
    f = lambda t, y: -y
    def err(fn, h): return abs(fn(f, 0, [1.0], h, 1.0)[1][-1, 0] - math.exp(-1))
    o2 = math.log(err(nm.rk2, 0.1) / err(nm.rk2, 0.05)) / math.log(2)
    o4 = math.log(err(nm.rk4, 0.1) / err(nm.rk4, 0.05)) / math.log(2)
    assert o2 == pytest.approx(2, abs=0.2) and o4 == pytest.approx(4, abs=0.3)


def test_interpolacion():
    xn = np.linspace(0, 2 * math.pi, 9); yn = np.sin(xn)
    T = nm.newton_table(xn, yn); S, _ = nm.spline_natural(xn, yn); xs = np.linspace(0, 2 * math.pi, 50)
    np.testing.assert_allclose(nm.lagrange(xn, yn, xs), nm.newton_eval(xn, T, xs), atol=1e-12)   # mismo polinomio
    assert np.max(np.abs(S(xs) - np.sin(xs))) < 0.05
    np.testing.assert_allclose(nm.lagrange(xn, yn, xn), yn, atol=1e-12)                          # pasa por los nodos
