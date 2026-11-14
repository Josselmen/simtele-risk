"""SIMTELE-RISK: simulador educativo de metodos numericos (Flask + SQLite).

Ejecutar:  python app.py   ->  http://127.0.0.1:5000
Todos los datos son sinteticos; el proyecto es academico. Entrega: Semana 1.
"""
import os

from flask import Flask, abort, render_template, request

from db import database
from modules import MODULES, RELACION, RELACION_HEADERS
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
        return dict(SCENARIOS=SCENARIOS, MODULES=MODULES, RELACION=RELACION, RELACION_HEADERS=RELACION_HEADERS)

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
        if not meta.get("implementado"):          # Semana 1: caso solo planteado (sin ejecucion numerica)
            return render_template("modulo.html", meta=meta, scenario=scenario, fields=fields,
                                   values={k: str(v) for k, v in defaults(meta["fields"], scenario).items()}, result=None, errors=[])
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
        return render_template("modulo.html", meta=meta, scenario=scenario, fields=fields, values=values, result=result, errors=errors)

    @app.route("/acerca")
    def acerca():
        return render_template("acerca.html")

    @app.errorhandler(404)
    def not_found(_e):
        return render_template("error.html", msg="La página solicitada no existe."), 404

    return app


if __name__ == "__main__":
    create_app().run(debug=True)
