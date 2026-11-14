"""SIMTELE-RISK: simulador educativo de metodos numericos (Flask + SQLite).

Ejecutar:  python app.py   ->  http://127.0.0.1:5000
Todos los datos son sinteticos; el proyecto es academico.
"""
import json
import os

from flask import Flask, abort, flash, redirect, render_template, request, url_for

from db import database
from modules import MODULES
from modules.common import SCENARIOS, ValidationError, clean, defaults, parse_params


def create_app(db_path=None):
    app = Flask(__name__)
    app.config["SECRET_KEY"] = os.environ.get("SECRET_KEY", "simtele-risk-dev")
    app.config["DB_PATH"] = db_path or database.DEFAULT_DB
    database.init_db(app.config["DB_PATH"])

    def fields_for(meta, scenario):
        return [f for f in meta["fields"] if "only" not in f or scenario in f["only"]]

    @app.context_processor
    def inject():
        return dict(SCENARIOS=SCENARIOS, MODULES=MODULES)

    @app.route("/")
    def index():
        return render_template("index.html")

    @app.route("/modulo/<int:mid>", methods=["GET", "POST"])
    def modulo(mid):
        mod = MODULES.get(mid) or abort(404)
        meta = mod.META
        scenario = request.values.get("escenario", "NORMAL")
        if scenario not in SCENARIOS:
            scenario = "NORMAL"
        fields = fields_for(meta, scenario)
        values = {k: str(v) for k, v in defaults(meta["fields"], scenario).items()}
        result, errors = None, []
        if request.method == "GET":
            for f in fields:                      # permite relanzar con ?v0=...&a0=...
                if f["name"] in request.args:
                    values[f["name"]] = request.args[f["name"]]
        else:
            values.update({f["name"]: request.form.get(f["name"], values.get(f["name"], "")) for f in fields})
            params, errors = parse_params(meta["fields"], request.form, scenario)
            if not errors:
                try:
                    result = clean(mod.run(scenario, params))
                except ValidationError as e:
                    errors.append(str(e))
                except (ValueError, ZeroDivisionError, OverflowError, FloatingPointError) as e:
                    errors.append(f"No se pudo calcular con esos parámetros: {e}")
            if result is not None and request.form.get("accion") == "guardar":
                sim_id = database.save_simulation(result, app.config["DB_PATH"])
                flash(f"Simulación #{sim_id} guardada en la base de datos.", "ok")
                return redirect(url_for("detalle", sim_id=sim_id))
        return render_template("modulo.html", meta=meta, scenario=scenario, fields=fields, values=values, result=result, errors=errors)

    @app.route("/historial")
    def historial():
        modulo_ = request.args.get("modulo") or None
        esc = request.args.get("escenario") or None
        texto = request.args.get("q") or None
        rows = database.list_simulations(app.config["DB_PATH"], modulo_, esc, texto)
        for r in rows:
            r["params"] = json.loads(r["parametros_json"])
            r["scenario"] = {"NORMAL": "NORMAL", "ESTRES": "ESTRES", "MAL_CONDICIONADO": "MAL"}[r["tipo"]]
        return render_template("historial.html", rows=rows, f_modulo=modulo_ or "", f_esc=esc or "", f_q=texto or "")

    @app.route("/historial/<int:sim_id>")
    def detalle(sim_id):
        sim = database.load_simulation(sim_id, app.config["DB_PATH"]) or abort(404)
        mod = MODULES[sim["caso_id"]]
        result, error = None, None
        try:                                       # se reconstruyen tablas y graficos desde los parametros guardados
            result = clean(mod.run(sim["scenario"], sim["params"]))
        except Exception as e:                     # noqa: BLE001
            error = str(e)
        stored = []
        by_name = {m["name"]: m for m in (result["methods"] if result else [])}
        for ej in sim["ejecuciones"]:
            res = ej["resumen"] or {}
            val = json.loads(res["valor_json"]) if res.get("valor_json") else None
            now = by_name.get(ej["metodo"])
            same = None
            if now is not None and isinstance(val, (int, float)) and isinstance(now["value"], (int, float)):
                same = abs(val - now["value"]) <= 1e-9 * max(1.0, abs(val))
            stored.append(dict(ej=ej, valor=val, ref=res.get("referencia"), err=json.loads(res["metrica_error_json"]) if res.get("metrica_error_json") else {}, same=same))
        query = "&".join(f"{k}={v}" for k, v in sim["params"].items())
        return render_template("detalle.html", sim=sim, result=result, error=error, stored=stored, meta=mod.META, query=query)

    @app.route("/historial/<int:sim_id>/eliminar", methods=["POST"])
    def eliminar(sim_id):
        database.delete_simulation(sim_id, app.config["DB_PATH"])
        flash(f"Simulación #{sim_id} eliminada.", "ok")
        return redirect(url_for("historial"))

    @app.route("/acerca")
    def acerca():
        return render_template("acerca.html")

    @app.errorhandler(404)
    def not_found(_e):
        return render_template("error.html", msg="La página solicitada no existe."), 404

    return app


if __name__ == "__main__":
    create_app().run(debug=True)
