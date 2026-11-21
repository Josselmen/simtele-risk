"""Modulo 4 - Reconstruccion de senales de sensores.

Senal sintetica y = sin t en [0, 2*pi]; se ocultan nodos (registros faltantes) y se
reconstruyen con Lagrange, Newton (diferencias divididas) y spline cubico natural.
Raices: cruce del umbral y = yu sobre el spline.
"""
import math
import time

import numpy as np

import numerics as nm
from .common import *  # noqa: F401,F403
from .common import plt, fig_b64, fmt, table, method, iters_from_hist, iter_table, ValidationError

MODULE_ID = 4
META = dict(
    id=4, slug="sensores", nombre="Señales de sensores", icono="📡",
    descripcion="¿Cómo reconstruir registros faltantes de un sensor y cuánto distorsiona el ruido la reconstrucción?",
    modelo="y(t) = sen t en [0, 2π], nodos equiespaciados; yᵢ = sen tᵢ + σ·ξᵢ (ξ ~ N(0,1)); se ocultan algunos nodos.",
    metodos="Lagrange · Newton (diferencias divididas) · Spline cúbico natural · Bisección · Newton · Secante",
    escenarios={
        "NORMAL": "11 nodos, 2 ocultos, sin ruido.",
        "ESTRES": "15 nodos, 4 ocultos y ruido σ = 0.05; se estudia el efecto del nivel de ruido.",
        "MAL": "25 nodos equiespaciados con ruido 1e-3: la constante de Lebesgue amplifica el ruido en la interpolación polinómica.",
    },
    limites=("Señal sintética y ruido gaussiano independiente: no representa sensores reales con saltos o ruido de otro tipo. "
             "No se extrapola fuera del intervalo muestreado ni se infiere la causa de la falla."),
    fields=[
        dict(name="n", label="Número de nodos n", unit="", type="int", min=7, max=40, defaults={"NORMAL": 11, "ESTRES": 15, "MAL": 25}),
        dict(name="hid", label="Índices ocultos (separados por coma)", unit="", type="text", defaults={"NORMAL": "3,7", "ESTRES": "2,5,8,11", "MAL": "6,12,18"}),
        dict(name="sigma", label="Ruido σ", unit="", type="float", min=0, max=0.5, defaults={"NORMAL": 0.0, "ESTRES": 0.05, "MAL": 0.001}),
        dict(name="seed", label="Semilla aleatoria", unit="", type="int", min=0, max=10 ** 6, defaults={"NORMAL": 42, "ESTRES": 42, "MAL": 11}),
        dict(name="trials", label="Realizaciones para estadísticas", unit="", type="int", min=1, max=200, defaults={"NORMAL": 30, "ESTRES": 50, "MAL": 30}),
        dict(name="yu", label="Umbral de cruce y_u", unit="", type="float", min=0.1, max=0.9, defaults={"NORMAL": 0.5, "ESTRES": 0.5, "MAL": 0.5}),
    ],
)
TW = 2 * math.pi
XS = np.linspace(0, TW, 400)


def _parse_hidden(txt, n):
    try:
        idx = sorted({int(s) for s in str(txt).replace(";", ",").split(",") if s.strip() != ""})
    except ValueError:
        raise ValidationError("Los índices ocultos deben ser enteros separados por coma (ejemplo: 3,7).")
    bad = [i for i in idx if i < 1 or i > n - 2]
    if bad:
        raise ValidationError(f"Los índices ocultos deben estar entre 1 y {n - 2} (los extremos no se ocultan). Inválidos: {bad}.")
    if n - len(idx) < 4:
        raise ValidationError("Deben quedar al menos 4 nodos disponibles para interpolar.")
    return idx


def _recon(tn, yn, hid):
    keep = [i for i in range(len(tn)) if i not in hid]
    xk, yk = tn[keep], yn[keep]
    T = nm.newton_table(xk, yk); S, _ = nm.spline_natural(xk, yk)
    return xk, yk, T, S, (lambda x: nm.lagrange(xk, yk, x)), (lambda x: nm.newton_eval(xk, T, x))


def run(scenario, p):
    n, sigma, seed, trials, yu = p["n"], p["sigma"], p["seed"], p["trials"], p["yu"]
    hid = _parse_hidden(p["hid"], n)
    rng = np.random.default_rng(seed)
    tn = np.linspace(0, TW, n); yclean = np.sin(tn); yn = yclean + sigma * rng.standard_normal(n)
    t0 = time.perf_counter(); xk, yk, Tn, S, L, N = _recon(tn, yn, hid); ms = 1e3 * (time.perf_counter() - t0)
    methods, tables = [], []
    th = tn[hid]; funcs = (("Lagrange", L), ("Newton", N), ("Spline cúbico", S))
    rows = []
    for t_ in th:
        row = [t_, math.sin(t_)]
        vals = [float(fn(t_)[0]) for _, fn in funcs]
        rows.append(row + vals + [abs(v - math.sin(t_)) for v in vals])
    tables.append(table("Reconstrucción en los puntos ocultos", ["t", "y real (sen t)", "Lagrange", "Newton", "Spline", "Error Lagrange", "Error Newton", "Error spline"], rows,
                        note=f"Nodos disponibles: {len(xk)} de {n}. Ruido σ = {fmt(sigma, 3)}."))
    dense = {nm_: float(np.max(np.abs(fn(XS) - np.sin(XS)))) for nm_, fn in funcs}
    dLN = float(np.max(np.abs(L(XS) - N(XS))))
    lam = nm.lebesgue(xk)
    for (nm_, fn), r in zip(funcs, rows):
        methods.append(method("Interpolación · " + nm_, float(fn(th[0])[0]), math.sin(th[0]), True, [], None, None, None, ms))
    tables.append(table("Error máximo sobre 400 puntos de [0, 2π]", ["Método", "Error máximo"], [[k, v] for k, v in dense.items()] + [["máx |Lagrange − Newton|", dLN], ["Constante de Lebesgue Λ (nodos disponibles)", lam]],
                        note="Lagrange y Newton construyen el mismo polinomio: solo difieren por redondeo."))
    Tt = nm.newton_table(xk, yk)
    rows_dd = [[i, fmt(xk[i], 5)] + [fmt(Tt[i, j], 5) if i + j < len(xk) else "" for j in range(5)] for i in range(len(xk))]
    tables.append(table("Tabla de diferencias divididas de Newton (primeras 5 columnas)", ["i", "tᵢ", "f[tᵢ]", "orden 1", "orden 2", "orden 3", "orden 4"], rows_dd))

    # ---- cruce de umbral sobre el spline
    g = lambda x: float(S(x)[0]) - yu
    dg = lambda x, e=1e-6: (g(x + e) - g(x - e)) / (2 * e)
    ref = math.asin(yu); TOL = 1e-10
    try:
        rt = {"Bisección": nm.bisection(g, 0, math.pi / 2, TOL, xref=ref), "Newton": nm.newton(g, dg, ref * 0.6, TOL, xref=ref), "Secante": nm.secant(g, ref * 0.6, ref * 1.2, TOL, xref=ref)}
        for k, r in rt.items():
            methods.append(method("Cruce y = y_u · " + k, r["root"], ref, r["converged"], iters_from_hist(r["iters"]), TOL, None, 100))
        tables.append(table(f"Instante de cruce del umbral y = {fmt(yu)} sobre el spline", ["Método", "Raíz", "Iteraciones", "Error absoluto"],
                            [[k, fmt(r["root"], 8), len(r["iters"]), fmt(r["iters"][-1]["err"], 3)] for k, r in rt.items()],
                            note=f"Referencia arcsen(y_u) = {fmt(ref, 8)}. El error restante proviene de la reconstrucción, no del método de raíces."))
    except ValueError as e:
        tables.append(table("Cruce de umbral", ["Resultado"], [[f"No se pudo calcular: {e}"]]))

    # ---- efecto del espaciamiento
    sp = []
    for nn in (7, 9, 11, 15, 21):
        t_ = np.linspace(0, TW, nn); y_ = np.sin(t_); hh = sorted({max(1, nn // 3), (2 * nn) // 3})
        _, _, _, S_, L_, N_ = _recon(t_, y_, hh)
        sp.append([nn, TW / (nn - 1)] + [float(np.max(np.abs(fn(XS) - np.sin(XS)))) for fn in (L_, N_, S_)])
    tables.append(table("Efecto del número de nodos (sin ruido, 2 nodos ocultos)", ["Nodos", "h = 2π/(n−1)", "Lagrange", "Newton", "Spline"], sp))

    # ---- efecto del ruido (promedio de realizaciones)
    nr = []
    rng2 = np.random.default_rng(seed + 1)
    for sg in sorted({0.0, 1e-3, 1e-2, 5e-2, 1e-1, sigma}):
        accL, accN, accS = [], [], []
        for _ in range(1 if sg == 0 else trials):
            yy = yclean + sg * rng2.standard_normal(n)
            _, _, _, S2, L2, N2 = _recon(tn, yy, hid)
            accL.append(np.max(np.abs(L2(XS) - np.sin(XS)))); accN.append(np.max(np.abs(N2(XS) - np.sin(XS)))); accS.append(np.max(np.abs(S2(XS) - np.sin(XS))))
        nr.append([sg, float(np.mean(accL)), float(np.mean(accN)), float(np.mean(accS)), (float(np.mean(accL)) / sg) if sg > 0 else None])
    tables.append(table(f"Error máximo medio según el nivel de ruido ({trials} realizaciones por nivel)", ["σ", "Lagrange", "Newton", "Spline", "Amplificación Lagrange (error/σ)"],
                        [[a, b, c_, d, (e if e is not None else "—")] for a, b, c_, d, e in nr],
                        note="Amplificación grande = el problema (nodos equiespaciados de grado alto) magnifica el ruido; el algoritmo no falla."))

    # ---- graficos
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
    ax[0].plot(XS, np.sin(XS), "k", lw=1, label="sen t")
    for (nm_, fn), ls in zip(funcs, ("-", "--", ":")):
        ax[0].plot(XS, fn(XS), ls, lw=1.2, label=nm_)
    ax[0].plot(xk, yk, "ko", ms=3, label="nodos"); ax[0].plot(th, np.sin(th), "rx", ms=7, label="ocultos")
    lo, hi_ = np.sin(XS).min() - 1.2, np.sin(XS).max() + 1.2; ax[0].set_ylim(lo, hi_)
    ax[0].set(xlabel="t", ylabel="y", title="Reconstrucción"); ax[0].legend(fontsize=6)
    for (nm_, fn), ls in zip(funcs, ("-", "--", ":")):
        ax[1].semilogy(XS, np.abs(fn(XS) - np.sin(XS)) + 1e-18, ls, lw=1, label=nm_)
    ax[1].set(xlabel="t", ylabel="|error|", title="Error punto a punto"); ax[1].legend(fontsize=7)
    pos = [r for r in nr if r[0] > 0]
    ax[2].loglog([r[0] for r in pos], [r[1] for r in pos], "o-", label="Lagrange/Newton"); ax[2].loglog([r[0] for r in pos], [r[3] for r in pos], "s-", label="Spline")
    ax[2].set(xlabel="σ", ylabel="error máx. medio", title="Error frente al ruido"); ax[2].legend(fontsize=7)
    charts = [dict(title="Reconstrucción, error y ruido", png=fig_b64(fig))]

    summary = table("Datos del experimento", ["Magnitud", "Valor"], [["Nodos / ocultos", f"{n} / {hid}"], ["Ruido σ", fmt(sigma, 3)], ["Constante de Lebesgue Λ", fmt(lam, 4)],
                                                                     ["Error máx. Lagrange", fmt(dense["Lagrange"], 3)], ["Error máx. spline", fmt(dense["Spline cúbico"], 3)]])
    amp_L = next((r[4] for r in nr if abs(r[0] - sigma) < 1e-15 and r[4]), None)
    num = (f"Lagrange y Newton coinciden (máx. diferencia {fmt(dLN, 2)}). Error máximo: Lagrange {fmt(dense['Lagrange'], 2)}, spline {fmt(dense['Spline cúbico'], 2)}. "
           f"La constante de Lebesgue de los nodos disponibles es Λ = {fmt(lam, 4)}." +
           (f" Con σ = {fmt(sigma, 3)} el polinomio amplifica el ruido ≈ {fmt(amp_L, 3)} veces." if amp_L else ""))
    prob = (f"Los {len(hid)} registros faltantes se reconstruyen con error ≈ {fmt(max(r[5] for r in rows), 2)} (Lagrange) o ≈ {fmt(max(r[7] for r in rows), 2)} (spline) sobre una señal de amplitud 1. " +
            ("Con ruido conviene el spline y evitar polinomios de grado alto en nodos equiespaciados; los extremos son los puntos menos confiables." if sigma > 0 else
             "Sin ruido el polinomio es más preciso; con ruido el orden suele invertirse."))
    return dict(module_id=MODULE_ID, scenario=scenario, params=p, summary=summary, tables=tables, charts=charts, methods=methods,
                interpretation=dict(numeric=num, problem=prob), limits=META["limites"])
