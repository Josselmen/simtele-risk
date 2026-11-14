"""Modulo 3 - Oscilaciones por viento.

Modelo: m x'' + c x' + k x = F0 sin(Om t)  ->  x' = v, v' = (F - c v - k x)/m.
Metodos: RK2 (punto medio) y RK4; trapecio/Simpson/Romberg para I = integral de x^2.
Escenario MAL: sensibilidad a k (problema) frente a crecimiento espurio de RK2 (algoritmo).
"""
import math
import time

import numpy as np

import numerics as nm
from .common import *  # noqa: F401,F403
from .common import plt, fig_b64, fmt, table, method, order_mean, ValidationError

MODULE_ID = 3
META = dict(
    id=3, slug="oscilaciones", nombre="Oscilaciones por viento", icono="🌬️",
    descripcion="¿Cómo cambia la oscilación de una cabina ante la fuerza del viento, el amortiguamiento y la rigidez?",
    modelo="m·x'' + c·x' + k·x = F₀·sen(Ω·t);  x' = v,  v' = (F − c·v − k·x)/m.  Exposición I = ∫ x(t)² dt.",
    metodos="RK2 (punto medio) · RK4 · Trapecio · Simpson · Romberg",
    escenarios={
        "NORMAL": "Vibración libre subamortiguada con solución analítica.",
        "ESTRES": "Viento armónico cerca de la resonancia, menor amortiguamiento, pasos grandes.",
        "MAL": "Casi sin amortiguamiento y T largo: sensibilidad a k y crecimiento espurio de RK2 (dos fenómenos distintos).",
    },
    limites=("Oscilador lineal de un grado de libertad con parámetros sintéticos: no representa el péndulo no lineal de una cabina real, "
             "viento turbulento de banda ancha, interacción con el cable ni amortiguamiento aerodinámico."),
    fields=[
        dict(name="m", label="Masa m", unit="kg", type="float", min=10, max=1e5, defaults={"NORMAL": 2000.0, "ESTRES": 2000.0, "MAL": 2000.0}),
        dict(name="k", label="Rigidez k", unit="N/m", type="float", min=10, max=1e6, defaults={"NORMAL": 8000.0, "ESTRES": 8000.0, "MAL": 8000.0}),
        dict(name="c", label="Amortiguamiento c", unit="N·s/m", type="float", min=0, max=1e5, defaults={"NORMAL": 400.0, "ESTRES": 200.0, "MAL": 16.0}),
        dict(name="F0", label="Amplitud del viento F₀", unit="N", type="float", min=0, max=1e5, defaults={"NORMAL": 0.0, "ESTRES": 600.0, "MAL": 0.0}),
        dict(name="Om", label="Frecuencia del viento Ω", unit="rad/s", type="float", min=0.01, max=50, defaults={"NORMAL": 1.8, "ESTRES": 1.8, "MAL": 1.8}),
        dict(name="x0", label="Posición inicial x₀", unit="m", type="float", min=-10, max=10, defaults={"NORMAL": 0.0, "ESTRES": 0.0, "MAL": 0.0}),
        dict(name="v0", label="Velocidad inicial v₀", unit="m/s", type="float", min=-10, max=10, defaults={"NORMAL": 1.0, "ESTRES": 0.0, "MAL": 1.0}),
        dict(name="T", label="Tiempo total T", unit="s", type="float", min=1, max=1000, defaults={"NORMAL": 20.0, "ESTRES": 60.0, "MAL": 300.0}),
        dict(name="kexp", label="Pasos N = 2^e (exponente e)", unit="", type="int", min=6, max=13, defaults={"NORMAL": 8, "ESTRES": 9, "MAL": 12}),
    ],
)
NREF = 2 ** 15


def _analytic(m, k, c, x0, v0):
    w0 = math.sqrt(k / m); z = c / (2 * math.sqrt(k * m))
    if z >= 1:
        return None
    wd = w0 * math.sqrt(1 - z * z)
    return lambda t: np.exp(-z * w0 * t) * (x0 * np.cos(wd * t) + (v0 + z * w0 * x0) / wd * np.sin(wd * t))


def run(scenario, p):
    m, k, c, F0, Om, x0, v0, T, kexp = (p[q] for q in ("m", "k", "c", "F0", "Om", "x0", "v0", "T", "kexp"))
    N = 2 ** int(kexp)
    w0 = math.sqrt(k / m); z = c / (2 * math.sqrt(k * m)); wd = w0 * math.sqrt(max(1 - z * z, 0.0))
    f = lambda t, y: np.array([y[1], (F0 * math.sin(Om * t) - c * y[1] - k * y[0]) / m])
    xa = _analytic(m, k, c, x0, v0) if F0 == 0 else None
    # ---- referencia
    t0 = time.perf_counter()
    if xa is not None:
        tref = np.linspace(0, T, NREF + 1); xref = xa(tref); ref_name = "solución analítica"
    else:
        tref, Yr = nm.rk4(f, 0, [x0, v0], T / NREF, T); xref = Yr[:, 0]; ref_name = f"RK4 de referencia con {NREF} pasos"
    href = T / NREF
    Iref = nm.simpson_arr(xref ** 2, href)
    amp_ref = float(np.max(np.abs(xref)))
    Xss = F0 / math.sqrt((k - m * Om ** 2) ** 2 + (c * Om) ** 2) if F0 > 0 else None
    methods, tables = [], []

    # ---- EDO con refinamiento
    Ns = [n for n in (N // 4, N // 2, N, 2 * N) if n >= 16 and n <= NREF]
    rows, e2, e4, hs = [], [], [], []
    traj = {}
    for n in Ns:
        h = T / n; hs.append(h); stride = NREF // n; row = [n, h]
        for name, fn in (("RK2", nm.rk2), ("RK4", nm.rk4)):
            t0 = time.perf_counter()
            with np.errstate(all="ignore"):
                Tn, Y = fn(f, 0, [x0, v0], h, T)
            ms = 1e3 * (time.perf_counter() - t0)
            err = float(np.max(np.abs(Y[:, 0] - xref[::stride]))); amp = float(np.max(np.abs(Y[:, 0])))
            with np.errstate(all="ignore"):
                I = nm.simpson_arr(Y[:, 0] ** 2, h)
            row += [err, amp, abs(I - Iref)]
            (e2 if name == "RK2" else e4).append(err); traj[(name, n)] = (Tn, Y)
            if n == N:
                methods.append(method(name, amp, amp_ref, True, step=h, ms=ms))
        rows.append(row)
    tables.append(table(f"EDO: RK2 y RK4 frente a la referencia ({ref_name})",
                        ["N", "h [s]", "RK2 error máx. x [m]", "RK2 amplitud máx. [m]", "RK2 error de I", "RK4 error máx. x [m]", "RK4 amplitud máx. [m]", "RK4 error de I"], rows,
                        note=f"Amplitud máxima de referencia: {fmt(amp_ref, 6)} m. Orden observado: RK2 ≈ {fmt(order_mean(hs[::-1], e2[::-1]) if False else order_mean(hs, e2), 3)}, RK4 ≈ {fmt(order_mean(hs, e4), 3)}. "
                             "El error de I usa Simpson sobre la trayectoria de cada método."))

    # ---- integracion de x^2
    irows, hh, et, es = [], [], [], []
    for kk in range(3, min(int(math.log2(NREF)), 11) + 1):
        n = 2 ** kk; stride = NREF // n; y = xref[::stride] ** 2; h = T / n
        T_ = nm.trap_arr(y, h); S_ = nm.simpson_arr(y, h); R_ = nm.romberg_arr(y, h)[-1, -1]
        irows.append([n, h, fmt(T_, 8), fmt(S_, 8), fmt(R_, 8), abs(T_ - Iref), abs(S_ - Iref), abs(R_ - Iref)])
        hh.append(h); et.append(abs(T_ - Iref)); es.append(abs(S_ - Iref))
    for name, val in (("Trapecio", T_), ("Simpson", S_), ("Romberg", R_)):
        methods.append(method("Integral de x² · " + name + f" (n={n})", float(val), Iref, True, step=h))
    tables.append(table(f"I = ∫ x(t)² dt sobre [0, {fmt(T, 4)}]", ["n", "h", "Trapecio", "Simpson", "Romberg", "Error trap.", "Error Simpson", "Error Romberg"], irows,
                        note=f"Referencia I = {fmt(Iref, 10)} m²·s."))

    # ---- sensibilidad y estabilidad (cuando hay referencia analitica)
    sens_amp = None
    if xa is not None:
        def E_rel(delta, Tm):
            kp = k * (1 + delta); xp = _analytic(m, kp, c, x0, v0)
            tt = np.linspace(0, Tm, 200001); a = xa(tt)
            return float(np.max(np.abs(xp(tt) - a)) / max(np.max(np.abs(a)), 1e-300))
        d0 = 1e-4
        srows = [[Tm, d0, E_rel(d0, Tm), E_rel(d0, Tm) / d0] for Tm in sorted({max(T / 30, 1), T / 6, T / 3, 2 * T / 3, T})]
        sens_amp = srows[-1][3]
        tables.append(table("Sensibilidad del problema: k → k(1+δ), δ = 1e-4", ["Horizonte T [s]", "δ", "Error relativo E = máx|Δx|/máx|x|", "Amplificación E/δ"], srows,
                            note="Es una propiedad del problema: no depende del método ni del paso."))
        drows = [[d, E_rel(d, T), E_rel(d, T) / d] for d in (1e-6, 1e-5, 1e-4, 1e-3)]
        tables.append(table(f"Linealidad de la sensibilidad (T = {fmt(T, 4)} s)", ["δ", "Error relativo E", "Amplificación E/δ"], drows))
        win = min(20.0, 0.2 * T); wref = float(np.max(np.abs(xa(np.linspace(T - win, T, 20001)))))
        grow = []
        for n in Ns:
            h = T / n; Rn = math.sqrt(1 + (w0 * h) ** 4 / 4)
            r2 = traj[("RK2", n)]; r4 = traj[("RK4", n)]
            w2 = r2[0] > T - win; w4 = r4[0] > T - win
            grow.append([h, Rn, float(np.max(np.abs(r2[1][w2, 0]))), float(np.max(np.abs(r4[1][w4, 0]))), wref])
        tables.append(table(f"Estabilidad numérica: amplitud máxima en los últimos {fmt(win, 3)} s", ["h [s]", "Factor |R| de RK2 por paso (sin amortiguar)", "RK2 [m]", "RK4 [m]", "Analítica [m]"], grow,
                            note="|R| > 1: RK2 amplifica la amplitud en cada paso; el efecto desaparece al reducir h (inestabilidad del algoritmo)."))

    # ---- graficos
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
    nplot = Ns[0]
    long_ = T > 100
    sel = (lambda a: a > T - 20) if long_ else (lambda a: a >= 0)
    tg = np.linspace(T - 20 if long_ else 0, T, 2000)
    ax[0].plot(tref[sel(tref)], xref[sel(tref)], "k", lw=1.2, label="referencia")
    for name, mk in (("RK2", "-"), ("RK4", "--")):
        Tn, Y = traj[(name, nplot)]; s_ = sel(Tn); ax[0].plot(Tn[s_], Y[s_, 0], mk, lw=1, label=f"{name} (N={nplot})")
    ax[0].set(xlabel="t [s]", ylabel="x [m]", title="Oscilación" + (" (últimos 20 s)" if long_ else "")); ax[0].legend(fontsize=7)
    ax[1].loglog(hs, [max(v, 1e-17) if np.isfinite(v) else np.nan for v in e2], "s--", label="RK2"); ax[1].loglog(hs, [max(v, 1e-17) if np.isfinite(v) else np.nan for v in e4], "o-", label="RK4")
    ax[1].set(xlabel="h [s]", ylabel="error máx. de x [m]", title="Error frente al paso"); ax[1].legend(fontsize=7)
    ax[2].loglog(hh, [max(v, 1e-18) for v in et], "^-", label="Trapecio"); ax[2].loglog(hh, [max(v, 1e-18) for v in es], "v-", label="Simpson")
    ax[2].set(xlabel="h [s]", ylabel="error de I", title="Integral de x²"); ax[2].legend(fontsize=7)
    charts = [dict(title="Trayectoria, error e integración", png=fig_b64(fig))]

    rows_s = [["Frecuencia natural ω₀ [rad/s]", fmt(w0, 5)], ["Razón de amortiguamiento ζ", fmt(z, 5)], ["Frecuencia amortiguada ω_d [rad/s]", fmt(wd, 6)],
              ["Amplitud máxima de la trayectoria [m]", fmt(amp_ref, 6)], ["I = ∫x² dt [m²·s]", fmt(Iref, 8)]]
    if Xss is not None:
        rows_s += [["Razón Ω/ω₀", fmt(Om / w0, 4)], ["Amplitud estacionaria teórica [m]", fmt(Xss, 5)], ["Deflexión estática F₀/k [m]", fmt(F0 / k, 4)], ["Amplificación dinámica", fmt(Xss / (F0 / k), 4)]]
    summary = table("Resultados de referencia", ["Magnitud", "Valor"], rows_s)
    num = (f"Con N = {N} (h = {fmt(T / N, 3)} s) el error máximo de x es {fmt(e4[Ns.index(N)], 2)} m con RK4 y {fmt(e2[Ns.index(N)], 2)} m con RK2; órdenes observados ≈ "
           f"{fmt(order_mean(hs, e2), 3)} y ≈ {fmt(order_mean(hs, e4), 3)}. ")
    if xa is not None:
        num += (f"La sensibilidad del problema a k es ≈ {fmt(sens_amp, 3)} (E/δ). Factor de RK2 por paso: {fmt(grow[Ns.index(N)][1], 8)}; amplitud calculada en la ventana final: "
                f"RK2 {fmt(grow[Ns.index(N)][2], 3)} m, RK4 {fmt(grow[Ns.index(N)][3], 3)} m, analítica {fmt(grow[Ns.index(N)][4], 3)} m.")
    else:
        num += "Con fuerza externa no hay sensibilidad analítica en esta versión; la referencia es RK4 con paso muy fino."
    prob = (f"Amplitud máxima {fmt(amp_ref, 4)} m; exposición acumulada I = {fmt(Iref, 4)} m²·s. " +
            (f"El viento está a {fmt(100 * Om / w0, 3)} % de la frecuencia natural y la amplitud estacionaria es {fmt(Xss / (F0 / k), 3)} veces la deflexión estática: "
             "cerca de la resonancia y con poco amortiguamiento la respuesta se amplifica mucho." if Xss is not None else
             f"La amplitud decae con constante de tiempo {fmt(1 / (z * w0), 3) if z > 0 else '∞'} s."))
    return dict(module_id=MODULE_ID, scenario=scenario, params=p, summary=summary, tables=tables, charts=charts, methods=methods,
                interpretation=dict(numeric=num, problem=prob), limits=META["limites"])
