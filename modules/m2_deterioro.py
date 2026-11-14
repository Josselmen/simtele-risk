"""Modulo 2 - Deterioro mecanico y mantenimiento.

Modelo: D(t) = D0 + al*t + be*t^2; t* con D(t*) = Du (referencia analitica) y
calibracion de coeficientes acoplados con un sistema lineal A x = b (Jacobi,
Gauss-Seidel, LU, gradiente conjugado). SEMANA 1: es el caso con los tres
escenarios ya probados; los metodos de raices quedan para la semana 2.
"""
import math
import time

import numpy as np

import numerics as nm
from .common import *  # noqa: F401,F403
from .common import plt, fig_b64, fmt, table, method, iters_from_lin, iter_table, ValidationError

MODULE_ID = 2
META = dict(
    id=2, slug="deterioro", nombre="Deterioro y mantenimiento", icono="🔧",
    descripcion="¿Cuándo alcanzaría un indicador hipotético de deterioro el umbral de inspección?",
    modelo="D(t) = D₀ + α·t + β·t²;  se resuelve D(t*) = D_u.  Calibración: A·x = b (SPD) con x* conocida.",
    metodos="Jacobi · Gauss-Seidel · LU · Gradiente conjugado",
    metodos_futuros="Bisección · Newton · Secante (D(t) − D_u = 0)",
    implementado=True,
    escenarios={
        "NORMAL": "Sistema SPD 4×4 diagonalmente dominante; deterioro nominal.",
        "ESTRES": "Mantenimiento retrasado (más desgaste, umbral menor) y sistema 20×20 con dominancia débil.",
        "MAL": "Matriz de Hilbert: Jacobi diverge, Gauss-Seidel es muy lento y las perturbaciones se amplifican.",
    },
    limites=("El indicador, el umbral y la parábola son hipotéticos: el deterioro real no es necesariamente cuadrático, "
             "depende del uso y es aleatorio. No se infieren plazos reales de mantenimiento."),
    fields=[
        dict(name="D0", label="Deterioro inicial D₀", unit="adim.", type="float", min=0, max=100, defaults={"NORMAL": 10.0, "ESTRES": 15.0, "MAL": 10.0}),
        dict(name="al", label="Tasa lineal α", unit="1/mes", type="float", min=0, max=10, defaults={"NORMAL": 0.8, "ESTRES": 1.1, "MAL": 0.8}),
        dict(name="be", label="Tasa cuadrática β", unit="1/mes²", type="float", min=0, max=1, defaults={"NORMAL": 0.02, "ESTRES": 0.035, "MAL": 0.02}),
        dict(name="Du", label="Umbral de inspección D_u", unit="adim.", type="float", min=1, max=1000, defaults={"NORMAL": 60.0, "ESTRES": 50.0, "MAL": 60.0}),
        dict(name="tplan", label="Inspección programada", unit="meses", type="float", min=1, max=600, defaults={"NORMAL": 24.0, "ESTRES": 24.0, "MAL": 24.0}),
        dict(name="n", label="Tamaño del sistema n", unit="", type="int", min=4, max=40, defaults={"NORMAL": 4, "ESTRES": 20, "MAL": 4}),
        dict(name="tol", label="Tolerancia", unit="", type="float", min=1e-14, max=1e-2, defaults={"NORMAL": 1e-10, "ESTRES": 1e-10, "MAL": 1e-10}),
    ],
)


def _system(scenario, n):
    if scenario == "NORMAL":
        A = nm.tridiag(n, 4.0, 1.0); xs = np.arange(1, n + 1, dtype=float)
        desc = f"tridiag(1, 4, 1) de tamaño {n}; x* = (1, 2, …, {n})"
    elif scenario == "ESTRES":
        A = nm.tridiag(n, 2.1, -1.0); xs = np.linspace(1, 2, n)
        desc = f"tridiag(−1, 2.1, −1) de tamaño {n}; x*ᵢ crece de 1 a 2"
    else:
        if n > 8:
            raise ValidationError("Para la matriz de Hilbert use n ≤ 8 (con n mayor la matriz es numéricamente singular).")
        A = nm.hilbert(n); xs = np.ones(n)
        desc = f"Hilbert {n}×{n}, Hᵢⱼ = 1/(i+j−1); x* = (1, …, 1)"
    return A, A @ xs, xs, desc


def run(scenario, p):
    D0, al, be, Du, tplan, n, tol = (p[k] for k in ("D0", "al", "be", "Du", "tplan", "n", "tol"))
    if D0 >= Du:
        raise ValidationError("D₀ debe ser menor que el umbral D_u.")
    if al <= 0 and be <= 0:
        raise ValidationError("Se necesita α > 0 o β > 0 para que el deterioro crezca.")
    methods, tables = [], []

    def tstar(D0_, al_, be_, Du_):
        if be_ < 1e-14:
            return (Du_ - D0_) / al_
        return (-al_ + math.sqrt(al_ * al_ - 4 * be_ * (D0_ - Du_))) / (2 * be_)
    ts = tstar(D0, al, be, Du)

    # ---- sensibilidad
    srows = []; names = {"al": "α", "be": "β", "D0": "D₀", "Du": "D_u"}; base = dict(D0=D0, al=al, be=be, Du=Du)
    for k in ("al", "be", "D0", "Du"):
        if base[k] == 0:
            continue
        for fac in (0.9, 1.1):
            q = dict(base); q[k] = base[k] * fac
            try:
                t2 = tstar(q["D0"], q["al"], q["be"], q["Du"])
            except ValueError:
                continue
            srows.append([names[k], f"{(fac - 1) * 100:+.0f} %", t2, t2 - ts, 100 * (t2 - ts) / ts, ((t2 - ts) / ts) / (fac - 1)])
    tables.append(table("Sensibilidad de t* a los parámetros (±10 %)", ["Parámetro", "Variación", "t* [meses]", "Δt* [meses]", "Δt*/t* [%]", "Coef. de sensibilidad S"], srows,
                        note="Referencia analítica de t*; la búsqueda numérica de t* (raíces) se implementa en la Semana 2."))

    # ---- sistema lineal de calibracion
    A, b, xs, desc = _system(scenario, n)
    sym = bool(np.allclose(A, A.T)); eigmin = float(np.linalg.eigvalsh(A).min()); spd = sym and eigmin > 0
    cond = float(np.linalg.cond(A)); dd = bool(np.all(2 * np.abs(np.diag(A)) > np.abs(A).sum(1)))
    D = np.diagflat(np.diag(A)); rhoJ = float(max(abs(np.linalg.eigvals(-np.linalg.inv(D) @ (A - D)))))
    tables.append(table("Sistema lineal de calibración", ["Propiedad", "Valor"], [
        ["Definición", desc], ["Simétrica", "Sí" if sym else "No"], ["Autovalor mínimo", fmt(eigmin, 4)], ["SPD (apta para CG)", "Sí" if spd else "No"],
        ["Diagonal dominante estricta", "Sí" if dd else "No"], ["ρ(T_Jacobi)", fmt(rhoJ, 4)], ["κ₂(A)", fmt(cond, 4)]]))
    lin = {}
    maxit = 20000
    t0 = time.perf_counter(); lin["Jacobi"] = nm.jacobi(A, b, tol, maxit, xs); ms_j = 1e3 * (time.perf_counter() - t0)
    t0 = time.perf_counter(); lin["Gauss-Seidel"] = nm.gauss_seidel(A, b, tol, maxit, xs); ms_g = 1e3 * (time.perf_counter() - t0)
    ms_c = None
    if spd:
        t0 = time.perf_counter(); lin["CG"] = nm.cg(A, b, tol, None, xs); ms_c = 1e3 * (time.perf_counter() - t0)
    t0 = time.perf_counter(); xlu = nm.lu_solve(A, b); ms_l = 1e3 * (time.perf_counter() - t0)
    rows = []
    for name, r in lin.items():
        it = r["iters"][-1]
        rows.append([name, r["n"], "Sí" if r["converged"] else f"No ({r['reason']})", fmt(it["res"], 3), fmt(it["err"], 3)])
        methods.append(method("Sistema lineal · " + name, float(it["x"][0]), float(xs[0]), r["converged"], iters_from_lin(r["iters"]), tol, None, maxit,
                              {"Jacobi": ms_j, "Gauss-Seidel": ms_g, "CG": ms_c}[name]))
    rows.append(["LU (directo)", "—", "Sí", fmt(nm.residual(A, b, xlu), 3), fmt(nm.err_abs(xlu, xs), 3)])
    methods.append(method("Sistema lineal · LU", float(xlu[0]), float(xs[0]), True, [], None, None, None, ms_l))
    if not spd:
        rows.append(["CG", "—", "No se ejecuta: la matriz no es SPD", "—", "—"])
    tables.append(table(f"Comparación de métodos (tolerancia relativa del residuo {tol:g}, x₀ = 0)", ["Método", "Iteraciones", "Convergió", "Residuo final", "Error ‖x − x*‖∞"], rows))
    bp = b.copy(); bp[0] *= (1 + 1e-6); xp = nm.lu_solve(A, bp)
    db_rel = abs(bp[0] - b[0]) / np.linalg.norm(b, np.inf); dx_rel = nm.err_abs(xp, xs) / np.linalg.norm(xs, np.inf)
    tables.append(table("Perturbación de b (b₁ cambia 1e-6 relativo)", ["Cambio relativo en x", "Amplificación", "Cota κ(A)"], [[dx_rel, dx_rel / db_rel, cond]]))
    for name in ("Jacobi", "Gauss-Seidel"):
        tables.append(iter_table(f"Iteraciones de {name}", iters_from_lin(lin[name]["iters"][:60]), 6, 0))

    # ---- graficos
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
    t = np.linspace(0, max(1.4 * ts, tplan * 1.2), 300)
    ax[0].plot(t, D0 + al * t + be * t * t, label="D(t)"); ax[0].axhline(Du, color="r", ls="--", label="umbral D_u")
    ax[0].axvline(tplan, color="gray", ls=":", label="inspección"); ax[0].plot([ts], [Du], "ko")
    ax[0].set(xlabel="t [meses]", ylabel="D", title="Deterioro y umbral (referencia analítica)"); ax[0].legend(fontsize=7)
    for name, r in lin.items():
        its = r["iters"]
        if len(its) > 1:
            ax[1].semilogy([i["k"] for i in its], [max(i["res"], 1e-18) for i in its], "-", lw=1, label=name)
    ax[1].set(xlabel="iteración k", ylabel="‖Ax − b‖∞", title=f"Residuo ({scenario}, κ={fmt(cond, 3)})"); ax[1].legend(fontsize=7)
    if scenario == "MAL":
        ax[1].set_xscale("log")
    charts = [dict(title="Deterioro y sistema lineal", png=fig_b64(fig))]

    summary = table("Resultados de referencia", ["Magnitud", "Valor", "Observación"], [
        ["Tiempo hasta el umbral t* [meses]", fmt(ts, 8), f"≈ {fmt(ts / 12, 3)} años (solución analítica)"],
        ["Raíz descartada (negativa) [meses]", fmt((-al - math.sqrt(al * al - 4 * be * (D0 - Du))) / (2 * be), 6) if be > 1e-14 else "—", "No tiene sentido físico"],
        ["Inspección programada [meses]", fmt(tplan), ""],
        ["Margen t* − t_plan [meses]", fmt(ts - tplan, 4), "Negativo: el umbral se alcanza antes de inspeccionar"],
    ])
    its_ = {k: v["n"] for k, v in lin.items()}
    gs = lin["Gauss-Seidel"]["iters"][-1]
    num = (f"En el sistema ({desc}), κ = {fmt(cond, 3)} y ρ(T_J) = {fmt(rhoJ, 3)}: " +
           ", ".join(f"{k} {'diverge' if not lin[k]['converged'] and lin[k]['reason'] == 'diverge' else str(v) + ' it.'}" for k, v in its_.items()) + ". " +
           (f"Gauss-Seidel termina con residuo {fmt(gs['res'], 2)} pero error {fmt(gs['err'], 2)} (error/residuo ≈ {fmt(gs['err'] / max(gs['res'], 1e-300), 3)}). " if scenario == "MAL" else "") +
           f"Una perturbación de 1e-6 en b se amplifica {fmt(dx_rel / db_rel, 3)} veces.")
    prob = (f"El indicador alcanza el umbral a los {fmt(ts, 4)} meses, " +
            (f"{fmt(ts - tplan, 3)} meses después de la inspección programada." if ts >= tplan else f"{fmt(tplan - ts, 3)} meses ANTES de la inspección programada.") +
            " " + (f"El parámetro más influyente es {max(srows, key=lambda r: abs(r[5]))[0]} (S ≈ {fmt(max(srows, key=lambda r: abs(r[5]))[5], 2)})." if srows else ""))
    return dict(module_id=MODULE_ID, scenario=scenario, params=p, summary=summary, tables=tables, charts=charts, methods=methods,
                interpretation=dict(numeric=num, problem=prob), limits=META["limites"])
