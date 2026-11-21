"""Modulo 1 - Detencion de emergencia de cabinas.

Modelo: x' = v, v' = -(a0 + c*v);  distancia total = v0*tr + x(t_det).
Metodos: RK2 y RK4 (EDO), trapecio/Simpson/Romberg (distancia = integral de v),
biseccion/Newton/secante (v(t) = 0). Escenario MAL CONDICIONADO: raiz doble.
"""
import math
import time

import numpy as np

import numerics as nm
from .common import *  # noqa: F401,F403
from .common import plt, fig_b64, fmt, table, method, iters_from_hist, iter_table, order_mean, ValidationError

MODULE_ID = 1
META = dict(
    id=1, slug="detencion", nombre="Detención de emergencia", icono="🚡",
    descripcion="¿Cuánto tiempo y distancia necesita una cabina para detenerse tras un frenado de emergencia?",
    modelo="x' = v,  v' = −(a₀ + c·v),  x(0)=0, v(0)=v₀.  Distancia total = v₀·tᵣ + x(t_det).",
    metodos="RK2 · RK4 · Trapecio · Simpson · Romberg · Bisección · Newton · Secante",
    escenarios={
        "NORMAL": "Frenado nominal con paso h adecuado.",
        "ESTRES": "Mayor velocidad, freno degradado, reacción lenta y pasos grandes.",
        "MAL": "v(t)=v₀(1−t/t_d)² tiene una raíz doble: la bisección no aplica y Newton pierde velocidad.",
    },
    limites=("Modelo unidimensional con parámetros sintéticos: no describe frenos reales, pendiente, viento, "
             "comodidad de pasajeros ni la dinámica del cable."),
    fields=[
        dict(name="v0", label="Velocidad inicial v₀", unit="m/s", type="float", min=0.1, max=30, defaults={"NORMAL": 5.0, "ESTRES": 6.0, "MAL": 5.0}),
        dict(name="a0", label="Desaceleración del freno a₀", unit="m/s²", type="float", min=0.05, max=10, defaults={"NORMAL": 1.0, "ESTRES": 0.6}, only=["NORMAL", "ESTRES"]),
        dict(name="c", label="Resistencia proporcional c", unit="1/s", type="float", min=0, max=5, defaults={"NORMAL": 0.2, "ESTRES": 0.2}, only=["NORMAL", "ESTRES"]),
        dict(name="tr", label="Tiempo de reacción tᵣ", unit="s", type="float", min=0, max=15, defaults={"NORMAL": 1.5, "ESTRES": 3.0}, only=["NORMAL", "ESTRES"]),
        dict(name="Dav", label="Distancia disponible", unit="m", type="float", min=1, max=500, defaults={"NORMAL": 20.0, "ESTRES": 20.0}, only=["NORMAL", "ESTRES"]),
        dict(name="h", label="Paso de RK (h)", unit="s", type="float", min=0.001, max=2, defaults={"NORMAL": 0.1, "ESTRES": 0.5}, only=["NORMAL", "ESTRES"]),
        dict(name="td", label="Tiempo de detención ideal t_d", unit="s", type="float", min=0.5, max=60, defaults={"MAL": 6.0}, only=["MAL"]),
    ],
)


def _analytic(v0, a0, c):
    if c < 1e-12:
        return (lambda t: v0 - a0 * t, lambda t: v0 * t - 0.5 * a0 * t * t, lambda t: -a0 + 0 * t, v0 / a0)
    K = v0 + a0 / c
    return (lambda t: K * np.exp(-c * t) - a0 / c,
            lambda t: K / c * (1 - np.exp(-c * t)) - a0 * t / c,
            lambda t: -K * c * np.exp(-c * t),
            math.log(1 + c * v0 / a0) / c)


def run(scenario, p):
    return _run_mal(p) if scenario == "MAL" else _run_normal(scenario, p)


def _run_normal(scenario, p):
    v0, a0, c, tr, Dav, h = (p[k] for k in ("v0", "a0", "c", "tr", "Dav", "h"))
    vf, xf, dvf, td = _analytic(v0, a0, c)
    d = float(xf(td)); d_react = v0 * tr; d_tot = d_react + d; margin = Dav - d_tot
    f = lambda t, y: np.array([y[1], -(a0 + c * y[1])])
    TOL = 1e-10
    methods, tables = [], []

    # ---- EDO con refinamiento de paso
    hs = [4 * h, 2 * h, h, h / 2]
    rows, e2s, e4s, ms = [], [], [], {}
    for hh in hs:
        row = [hh]
        for name, fn in (("RK2", nm.rk2), ("RK4", nm.rk4)):
            t0 = time.perf_counter(); T, Y = fn(f, 0, [0, v0], hh, td); ms[(name, hh)] = 1e3 * (time.perf_counter() - t0)
            row += [abs(Y[-1, 0] - d), abs(Y[-1, 1])]
            (e2s if name == "RK2" else e4s).append(abs(Y[-1, 0] - d))
        rows.append(row)
    for name in ("RK2", "RK4"):
        fn = nm.rk2 if name == "RK2" else nm.rk4
        T, Y = fn(f, 0, [0, v0], h, td)
        methods.append(method(name, float(Y[-1, 0]), d, True, step=h, ms=ms[(name, h)]))
    tables.append(table("EDO: RK2 (punto medio) y RK4 integrados hasta t_det",
                        ["h [s]", "RK2 error x(t_det) [m]", "RK2 |v(t_det)| [m/s]", "RK4 error x(t_det) [m]", "RK4 |v(t_det)| [m/s]"], rows,
                        note=f"Orden observado: RK2 ≈ {fmt(order_mean(hs, e2s), 3)}, RK4 ≈ {fmt(order_mean(hs, e4s), 3)} (teóricos 2 y 4)."))

    # ---- raices de v(t)=0
    hi = max(10.0, 3 * td)
    rt = {}
    for name, call in (("Bisección", lambda: nm.bisection(vf, 0, hi, TOL, xref=td)),
                       ("Newton", lambda: nm.newton(vf, dvf, 0.0, TOL, xref=td)),
                       ("Secante", lambda: nm.secant(vf, 0.0, 1.0, TOL, xref=td))):
        t0 = time.perf_counter(); r = call(); ms_ = 1e3 * (time.perf_counter() - t0)
        rt[name] = r
        methods.append(method("Raíz v(t)=0 · " + name, r["root"], td, r["converged"], iters_from_hist(r["iters"]), TOL, None, 100, ms_))
    tables.append(table("Raíces de v(t) = 0 (instante de detención)", ["Método", "Raíz [s]", "Iteraciones", "Convergió", "Error absoluto", "|v(raíz)|"],
                        [[k, fmt(v["root"], 10), len(v["iters"]), "Sí" if v["converged"] else "No", fmt(v["iters"][-1]["err"], 3), fmt(abs(v["iters"][-1]["fx"]), 3)] for k, v in rt.items()],
                        note=f"Referencia analítica: t_det = {fmt(td, 8)} s. Tolerancia {TOL:g}."))
    tables.append(iter_table("Iteraciones de Newton", rt["Newton"]["iters"]))
    tables.append(iter_table("Iteraciones de la bisección", rt["Bisección"]["iters"], 6, 2))

    # ---- integracion de v(t)
    irows, ns_h, et, es = [], [], [], []
    for k in range(1, 7):
        n = 2 ** k; x, y, hh = nm.grid(vf, 0, td, n)
        T_ = nm.trap_arr(y, hh); S_ = nm.simpson_arr(y, hh); R_ = nm.romberg_arr(y, hh)[-1, -1]
        irows.append([n, hh, fmt(T_, 8), fmt(S_, 8), fmt(R_, 8), abs(T_ - d), abs(S_ - d), abs(R_ - d)])
        ns_h.append(hh); et.append(abs(T_ - d)); es.append(abs(S_ - d))
    for name, val in (("Trapecio", T_), ("Simpson", S_), ("Romberg", R_)):
        methods.append(method("Integral de v(t) · " + name + " (n=64)", float(val), d, True, step=ns_h[-1]))
    tables.append(table("Distancia de frenado = ∫ v(t) dt", ["n", "h", "Trapecio", "Simpson", "Romberg", "Error trap.", "Error Simpson", "Error Romberg"], irows,
                        note=f"Orden observado: trapecio ≈ {fmt(order_mean(ns_h, et), 3)}, Simpson ≈ {fmt(order_mean(ns_h[:4], es[:4]), 3)}. Referencia: {fmt(d, 8)} m."))

    # ---- graficos
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
    t = np.linspace(0, td, 200)
    ax[0].plot(t, vf(t), "k", label="v(t) analítica")
    for name, mk in (("RK2", "s"), ("RK4", "o")):
        fn = nm.rk2 if name == "RK2" else nm.rk4
        T, Y = fn(f, 0, [0, v0], 4 * h, td); ax[0].plot(T, Y[:, 1], mk, ms=4, mfc="none", label=f"{name} (h={fmt(4 * h, 3)} s)")
    ax2 = ax[0].twinx(); ax2.plot(t, xf(t), "tab:red", lw=1); ax2.set_ylabel("x(t) [m]", color="tab:red"); ax2.grid(False)
    ax[0].set(xlabel="t [s]", ylabel="v [m/s]", title="Frenado: velocidad (negro) y distancia (rojo)"); ax[0].legend(fontsize=7)
    ax[1].loglog(hs, [max(v, 1e-17) for v in e2s], "s--", label="RK2: error en x(t_det)")
    ax[1].loglog(hs, [max(v, 1e-17) for v in e4s], "o-", label="RK4: error en x(t_det)")
    ax[1].loglog([td / r[0] for r in irows], [max(r[5], 1e-17) for r in irows], "^:", label="Trapecio (vs h)")
    ax[1].loglog([td / r[0] for r in irows], [max(r[6], 1e-17) for r in irows], "v:", label="Simpson (vs h)")
    ax[1].set(xlabel="paso h [s]", ylabel="error absoluto [m]", title="Error frente al paso"); ax[1].legend(fontsize=7)
    charts = [dict(title="Solución y error", png=fig_b64(fig))]

    summary = table("Resultados de referencia (analíticos)", ["Magnitud", "Valor", "Observación"], [
        ["Tiempo de detención t_det [s]", fmt(td, 8), f"Con c = 0 sería {fmt(v0 / a0, 4)} s"],
        ["Distancia de frenado [m]", fmt(d, 8), f"Con c = 0 sería {fmt(v0 ** 2 / (2 * a0), 4)} m"],
        ["Distancia en reacción v₀·tᵣ [m]", fmt(d_react, 5), "Velocidad constante antes de frenar"],
        ["Distancia total [m]", fmt(d_tot, 6), ""],
        ["Margen frente a la distancia disponible [m]", fmt(margin, 4), "Positivo: alcanza · Negativo: se excede"],
    ])
    num = (f"Con h = {fmt(h, 3)} s, RK4 comete {fmt(e4s[2], 2)} m de error en x(t_det) frente a {fmt(e2s[2], 2)} m de RK2; los órdenes observados son "
           f"≈ {fmt(order_mean(hs, e2s), 3)} y ≈ {fmt(order_mean(hs, e4s), 3)}. Para v(t)=0, Newton usó {len(rt['Newton']['iters'])} iteraciones, la secante {len(rt['Secante']['iters'])} "
           f"y la bisección {len(rt['Bisección']['iters'])}. Con 64 subintervalos los errores de integración son {fmt(et[-1], 2)} (trapecio) y {fmt(es[-1], 2)} (Simpson).")
    prob = (f"La cabina se detiene en {fmt(td, 4)} s y recorre {fmt(d, 4)} m al frenar. Con {fmt(d_react, 3)} m de reacción la distancia total es {fmt(d_tot, 4)} m, "
            + (f"que queda dentro de los {fmt(Dav, 4)} m disponibles (margen {fmt(margin, 3)} m)." if margin >= 0 else
               f"que EXCEDE en {fmt(-margin, 3)} m los {fmt(Dav, 4)} m disponibles."))
    if margin < 0:
        prob += f" Ese déficit es {fmt(-margin / max(e2s[0], 1e-12), 3)} veces mayor que el mayor error numérico de la tabla, así que la conclusión no depende del método."
    return dict(module_id=MODULE_ID, scenario=scenario, params=p, summary=summary, tables=tables, charts=charts, methods=methods,
                interpretation=dict(numeric=num, problem=prob), limits=META["limites"])


def _run_mal(p):
    v0, td = p["v0"], p["td"]
    vm = lambda t: v0 * (1 - t / td) ** 2
    dvm = lambda t: -2 * v0 * (1 - t / td) / td
    TOL = 1e-10; methods, tables = [], []
    try:
        nm.bisection(vm, 0, 1.5 * td); bis = "aplicable"
    except ValueError as e:
        bis = f"NO aplicable: {e}"
    newton = nm.newton(vm, dvm, 0.2 * td, TOL, xref=td, use_res=False)
    secant = nm.secant(vm, 0.0, 0.2 * td, TOL, xref=td, use_res=False)
    mod = nm.newton(vm, dvm, 0.2 * td, TOL, xref=td, m=2.0, use_res=True)
    resid = nm.newton(vm, dvm, 0.2 * td, TOL, xref=td, use_res=True)
    for name, r in (("Newton (m=1)", newton), ("Secante", secant), ("Newton modificado (m=2)", mod)):
        methods.append(method("Raíz doble · " + name, r["root"], td, r["converged"], iters_from_hist(r["iters"]), TOL, None, 100))
    tables.append(table("Búsqueda de la raíz doble t = t_d", ["Método", "Iteraciones", "Convergió", "Error |t − t_d| [s]", "|v(t)| [m/s]"],
                        [["Bisección", "—", "—", bis, "—"]] +
                        [[n_, len(r["iters"]), "Sí" if r["converged"] else "No", fmt(r["iters"][-1]["err"], 3), fmt(abs(r["iters"][-1]["fx"]), 3)]
                         for n_, r in (("Newton (m=1)", newton), ("Secante", secant), ("Newton modificado (m=2)", mod))],
                        note="Criterio de parada: tamaño del paso < 1e-10 (en el Newton modificado, residuo)."))
    nn = newton["iters"]; ratios = [nn[i + 1]["err"] / nn[i]["err"] for i in range(min(5, len(nn) - 1)) if nn[i]["err"] > 0]
    ratio = sum(ratios) / len(ratios) if ratios else float("nan")
    tables.append(iter_table("Iteraciones de Newton (convergencia lineal)", newton["iters"], 8, 2))
    tables.append(table("Parada por residuo: el residuo pequeño no garantiza error pequeño",
                        ["Criterio", "Iteraciones", "Error real |t − t_d| [s]", "|v(t)| al parar"],
                        [["|v(t)| < 1e-10", len(resid["iters"]), fmt(resid["iters"][-1]["err"], 3), fmt(abs(resid["iters"][-1]["fx"]), 3)],
                         ["paso < 1e-10", len(newton["iters"]), fmt(newton["iters"][-1]["err"], 3), fmt(abs(newton["iters"][-1]["fx"]), 3)]]))
    srows = []
    for eps in (1e-12, 1e-10, 1e-8, 1e-6):
        b = nm.bisection(lambda t: vm(t) - eps, 0, td, 1e-13); shift = td - b["root"]
        srows.append([eps, shift, td * math.sqrt(eps / v0), shift / eps])
    tables.append(table("Sensibilidad: se resuelve v(t) − ε = 0", ["ε [m/s]", "Desplazamiento de la raíz [s]", "Predicción t_d·√(ε/v₀) [s]", "Amplificación (desplaz./ε)"], srows,
                        note="El desplazamiento crece como √ε, no como ε: el problema está mal condicionado."))
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
    t = np.linspace(0, 1.4 * td, 300); ax[0].plot(t, vm(t), "tab:green"); ax[0].plot([td], [0], "ko")
    ax[0].set(xlabel="t [s]", ylabel="v [m/s]", title="v(t) = v₀(1 − t/t_d)²  (raíz doble en t_d)")
    for name, r in (("Newton (m=1)", newton), ("Secante", secant), ("Newton modificado", mod)):
        ax[1].semilogy([i["k"] for i in r["iters"]], [max(i["err"], 1e-17) for i in r["iters"]], "o-", ms=3, label=name)
    ax[1].set(xlabel="iteración k", ylabel="|t_k − t_d| [s]", title="Convergencia hacia la raíz doble"); ax[1].legend(fontsize=7)
    summary = table("Datos del escenario", ["Magnitud", "Valor"], [["v₀ [m/s]", fmt(v0)], ["t_d [s]", fmt(td)], ["Distancia analítica v₀·t_d/3 [m]", fmt(v0 * td / 3, 6)], ["Bisección", bis]])
    num = (f"Bisección: {bis}. Newton converge, pero de forma lineal con razón ≈ {fmt(ratio, 3)}: necesita {len(newton['iters'])} iteraciones. "
           f"La secante usa {len(secant['iters'])} y Newton modificado con la multiplicidad correcta converge en {len(mod['iters'])}. "
           f"Si se detiene por residuo, el error real es {fmt(resid['iters'][-1]['err'], 3)} s. Una perturbación de 1e-12 m/s desplaza la raíz {fmt(srows[0][1], 3)} s.")
    prob = ("Un frenado que suaviza la desaceleración hasta cero al detenerse produce una raíz doble: el instante de detención es muy sensible a pequeñas incertidumbres de la velocidad medida. "
            "Esto es mal condicionamiento del problema, no un fallo de convergencia de Newton.")
    return dict(module_id=MODULE_ID, scenario="MAL", params=p, summary=summary, tables=tables, charts=[dict(title="Raíz doble", png=fig_b64(fig))], methods=methods,
                interpretation=dict(numeric=num, problem=prob), limits=META["limites"])
