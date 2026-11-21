"""Utilidades comunes de los modulos: formato, validacion de parametros, tablas y graficos."""
import base64
import io
import math

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SCENARIOS = {"NORMAL": "NORMAL", "ESTRES": "ESTRÉS", "MAL": "MAL CONDICIONADO"}


class ValidationError(Exception):
    pass


def fmt(x, n=4):
    """Formato general: fijo o cientifico segun la magnitud."""
    if isinstance(x, (bool, np.bool_)):
        return "Sí" if x else "No"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if x is None:
        return "—"
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    if x != x:
        return "n/a"
    if x == 0:
        return "0"
    if abs(x) < 1e-3 or abs(x) >= 1e5:
        return f"{x:.{n - 1}e}"
    return f"{x:.{n}g}"


def table(title, headers, rows, note=""):
    return {"title": title, "headers": headers, "rows": [[fmt(c) if not isinstance(c, str) else c for c in r] for r in rows], "note": note}


def clean(o):
    """Convierte tipos NumPy a tipos nativos para JSON/plantillas."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def _dflt(f, scenario):
    d = f["defaults"]
    return d[scenario] if scenario in d else next(iter(d.values()))


def parse_params(fields, raw, scenario):
    """Valida los parametros recibidos del formulario. Devuelve (params, errores)."""
    params, errors = {}, []
    for f in fields:
        if "only" in f and scenario not in f["only"]:
            continue
        name = f["name"]
        txt = raw.get(name)
        if txt is None or str(txt).strip() == "":
            txt = _dflt(f, scenario)
        try:
            if f["type"] == "text":
                val = str(txt).strip()
            elif f["type"] == "int":
                val = int(float(str(txt).replace(",", ".")))
            else:
                val = float(str(txt).replace(",", "."))
        except ValueError:
            errors.append(f"«{f['label']}» debe ser un número válido (recibido: {txt!r}).")
            continue
        if f["type"] != "text":
            if not math.isfinite(val):
                errors.append(f"«{f['label']}» debe ser finito.")
                continue
            if "min" in f and val < f["min"]:
                errors.append(f"«{f['label']}» debe ser ≥ {f['min']} (recibido: {val:g}).")
            if "max" in f and val > f["max"]:
                errors.append(f"«{f['label']}» debe ser ≤ {f['max']} (recibido: {val:g}).")
            if f.get("positive") and val <= 0:
                errors.append(f"«{f['label']}» debe ser mayor que 0.")
        params[name] = val
    return params, errors


def defaults(fields, scenario):
    out = {}
    for f in fields:
        if "only" in f and scenario not in f["only"]:
            continue
        out[f["name"]] = _dflt(f, scenario)
    return out


def fig_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def new_fig(w=6.4, h=3.4):
    plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.3})
    return plt.subplots(figsize=(w, h))


def method(name, value, ref=None, converged=True, iters=None, tol=None, step=None, max_iter=None, ms=None):
    """Registro uniforme de un metodo ejecutado (lo que se guarda en la base de datos)."""
    m = dict(name=name, value=value, ref=ref, converged=bool(converged), iters=iters or [], tol=tol, step=step, max_iter=max_iter, ms=ms)
    try:
        if ref is not None and value is not None and not isinstance(value, (list, str)):
            m["err_abs"] = abs(float(value) - float(ref))
            m["err_rel"] = m["err_abs"] / abs(float(ref)) if abs(float(ref)) > 1e-14 else None
    except (TypeError, ValueError):
        pass
    return m


def iters_from_hist(hist, key="x"):
    """Convierte el historial de numerics.roots a la estructura guardable."""
    out = []
    for h in hist:
        out.append(dict(k=h["k"], x=h[key], res=abs(h.get("fx", 0.0)), err=h.get("err"), err_rel=None))
    return out


def iters_from_lin(hist):
    return [dict(k=h["k"], x=[float(v) for v in h["x"]], res=h["res"], err=h["err"], err_rel=None) for h in hist]


def iter_table(title, hist, first=8, last=2):
    n = len(hist)
    idx = list(range(min(first, n)))
    if n > first:
        idx += list(range(max(first, n - last), n))
    rows = []
    for i in idx:
        h = hist[i]
        x = h["x"]
        xs = fmt(x, 10) if not isinstance(x, (list, tuple)) else "(" + ", ".join(fmt(v, 6) for v in x[:6]) + (", …)" if len(x) > 6 else ")")
        rows.append([h["k"], xs, fmt(h.get("res"), 3), fmt(h.get("err"), 3)])
    return table(title, ["k", "Aproximación", "Residuo", "Error"], rows,
                 note="Se muestran las primeras y las últimas iteraciones." if n > first + last else "")


def order_mean(hs, es):
    o = [math.log(es[i] / es[i + 1]) / math.log(hs[i] / hs[i + 1]) for i in range(min(len(hs), len(es)) - 1) if es[i] > 0 and es[i + 1] > 0]
    return sum(o) / len(o) if o else float("nan")
