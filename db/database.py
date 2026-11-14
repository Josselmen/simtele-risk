"""Capa de datos SQLite de SIMTELE-RISK: crear, guardar, listar y recuperar simulaciones."""
import json
import os
import sqlite3

BASE = os.path.dirname(os.path.abspath(__file__))
SCHEMA_PATH = os.path.join(BASE, "schema.sql")
DEFAULT_DB = os.path.join(BASE, "simtele.db")

CASOS = [
    (1, "Detencion de emergencia"),
    (2, "Deterioro mecanico y mantenimiento"),
    (3, "Oscilaciones por viento"),
    (4, "Reconstruccion de senales de sensores"),
]


def connect(path=None):
    con = sqlite3.connect(path or DEFAULT_DB)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys = ON")
    return con


def init_db(path=None):
    """Crea las tablas si no existen y carga los 4 casos de estudio."""
    path = path or DEFAULT_DB
    nuevo = not os.path.exists(path)
    con = connect(path)
    existe = con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='caso_estudio'").fetchone()
    if not existe:
        with open(SCHEMA_PATH, encoding="utf-8") as fh:
            con.executescript(fh.read())
    for cid, nombre in CASOS:
        con.execute("INSERT OR IGNORE INTO caso_estudio(id, nombre, modelo_version) VALUES (?,?,?)", (cid, nombre, "v1.0"))
    con.commit()
    con.close()
    return nuevo


def _num(x):
    try:
        x = float(x)
        return x if x == x and abs(x) != float("inf") else None
    except (TypeError, ValueError):
        return None


def save_simulation(result, path=None, usuario="equipo"):
    """Guarda una simulacion completa. `result` es lo que devuelve modules.<m>.run()."""
    con = connect(path)
    cur = con.cursor()
    tipo = {"NORMAL": "NORMAL", "ESTRES": "ESTRES", "MAL": "MAL_CONDICIONADO"}[result["scenario"]]
    cur.execute("INSERT INTO escenario(caso_id,tipo,parametros_json,fuente_datos) VALUES (?,?,?,?)",
                (result["module_id"], tipo, json.dumps(result["params"]), "sintetico"))
    esc = cur.lastrowid
    cur.execute("INSERT INTO simulacion(escenario_id, usuario, estado) VALUES (?,?,?)", (esc, usuario, "completada"))
    sim = cur.lastrowid
    for m in result["methods"]:
        cur.execute("INSERT INTO ejecucion_metodo(simulacion_id,metodo,tolerancia,paso,max_iter,convergio,tiempo_ms) VALUES (?,?,?,?,?,?,?)",
                    (sim, m["name"], _num(m.get("tol")), _num(m.get("step")), m.get("max_iter"), int(bool(m.get("converged"))), _num(m.get("ms"))))
        ej = cur.lastrowid
        for it in m.get("iters", []):
            cur.execute("INSERT INTO iteracion_resultado(ejecucion_id,numero,aproximacion_json,residual,error_abs,error_rel) VALUES (?,?,?,?,?,?)",
                        (ej, it["k"], json.dumps(it["x"]), _num(it.get("res")), _num(it.get("err")), _num(it.get("err_rel"))))
        cur.execute("INSERT INTO resultado_resumen(ejecucion_id,referencia,valor_json,metrica_error_json,interpretacion) VALUES (?,?,?,?,?)",
                    (ej, _num(m.get("ref")), json.dumps(m.get("value")),
                     json.dumps({"err_abs": _num(m.get("err_abs")), "err_rel": _num(m.get("err_rel")), "iteraciones": len(m.get("iters", []))}),
                     result.get("interpretation", {}).get("problem", "")[:500]))
    con.commit()
    con.close()
    return sim


def list_simulations(path=None, modulo=None, escenario=None, texto=None):
    con = connect(path)
    q = """SELECT s.id, s.fecha, s.usuario, c.id AS caso_id, c.nombre AS caso, e.tipo, e.parametros_json,
                  (SELECT COUNT(*) FROM ejecucion_metodo x WHERE x.simulacion_id = s.id) AS n_metodos
           FROM simulacion s JOIN escenario e ON e.id = s.escenario_id JOIN caso_estudio c ON c.id = e.caso_id WHERE 1=1"""
    args = []
    if modulo:
        q += " AND c.id = ?"; args.append(int(modulo))
    if escenario:
        q += " AND e.tipo = ?"; args.append(escenario)
    if texto:
        q += " AND (c.nombre LIKE ? OR e.parametros_json LIKE ? OR s.usuario LIKE ?)"; args += [f"%{texto}%"] * 3
    q += " ORDER BY s.id DESC"
    rows = [dict(r) for r in con.execute(q, args).fetchall()]
    con.close()
    return rows


def load_simulation(sim_id, path=None):
    """Recupera parametros, ejecuciones, iteraciones y resumenes guardados."""
    con = connect(path)
    r = con.execute("""SELECT s.id, s.fecha, s.usuario, s.estado, c.id AS caso_id, c.nombre AS caso, e.tipo, e.parametros_json
                       FROM simulacion s JOIN escenario e ON e.id = s.escenario_id JOIN caso_estudio c ON c.id = e.caso_id
                       WHERE s.id = ?""", (sim_id,)).fetchone()
    if r is None:
        con.close()
        return None
    sim = dict(r)
    sim["params"] = json.loads(sim.pop("parametros_json"))
    sim["scenario"] = {"NORMAL": "NORMAL", "ESTRES": "ESTRES", "MAL_CONDICIONADO": "MAL"}[sim["tipo"]]
    sim["ejecuciones"] = []
    for ej in con.execute("SELECT * FROM ejecucion_metodo WHERE simulacion_id = ? ORDER BY id", (sim_id,)).fetchall():
        d = dict(ej)
        d["iteraciones"] = [dict(i) for i in con.execute(
            "SELECT numero, aproximacion_json, residual, error_abs, error_rel FROM iteracion_resultado WHERE ejecucion_id = ? ORDER BY numero", (ej["id"],)).fetchall()]
        rs = con.execute("SELECT referencia, valor_json, metrica_error_json, interpretacion FROM resultado_resumen WHERE ejecucion_id = ?", (ej["id"],)).fetchone()
        d["resumen"] = dict(rs) if rs else None
        sim["ejecuciones"].append(d)
    con.close()
    return sim


def delete_simulation(sim_id, path=None):
    con = connect(path)
    ejs = [r[0] for r in con.execute("SELECT id FROM ejecucion_metodo WHERE simulacion_id = ?", (sim_id,)).fetchall()]
    for ej in ejs:
        con.execute("DELETE FROM iteracion_resultado WHERE ejecucion_id = ?", (ej,))
        con.execute("DELETE FROM resultado_resumen WHERE ejecucion_id = ?", (ej,))
    con.execute("DELETE FROM ejecucion_metodo WHERE simulacion_id = ?", (sim_id,))
    esc = con.execute("SELECT escenario_id FROM simulacion WHERE id = ?", (sim_id,)).fetchone()
    con.execute("DELETE FROM simulacion WHERE id = ?", (sim_id,))
    if esc:
        con.execute("DELETE FROM escenario WHERE id = ? AND NOT EXISTS (SELECT 1 FROM simulacion WHERE escenario_id = ?)", (esc[0], esc[0]))
    con.commit()
    con.close()
