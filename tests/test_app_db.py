"""Pruebas de los modulos, la base de datos y las rutas web."""
import json

import pytest

from app import create_app
from db import database
from modules import MODULES
from modules.common import ValidationError, defaults, parse_params


@pytest.fixture()
def client(tmp_path):
    app = create_app(str(tmp_path / "test.db"))
    app.config["TESTING"] = True
    return app.test_client(), app


@pytest.mark.parametrize("mid", [1, 2, 3, 4])
@pytest.mark.parametrize("sc", ["NORMAL", "ESTRES", "MAL"])
def test_modulo_escenario(mid, sc):
    mod = MODULES[mid]
    r = mod.run(sc, defaults(mod.META["fields"], sc))
    assert r["tables"] and r["charts"] and r["methods"]
    assert r["interpretation"]["numeric"] and r["interpretation"]["problem"]


def test_modulo1_valores_analiticos():
    r = MODULES[1].run("NORMAL", defaults(MODULES[1].META["fields"], "NORMAL"))
    assert "3.4657" in r["summary"]["rows"][0][1] and "7.6713" in r["summary"]["rows"][1][1]


def test_validacion_de_parametros():
    _, errs = parse_params(MODULES[1].META["fields"], {"v0": "-3", "h": "abc"}, "NORMAL")
    assert len(errs) == 2
    with pytest.raises(ValidationError):
        MODULES[2].run("NORMAL", {**defaults(MODULES[2].META["fields"], "NORMAL"), "D0": 80.0})   # D0 >= Du
    with pytest.raises(ValidationError):
        MODULES[4].run("NORMAL", {**defaults(MODULES[4].META["fields"], "NORMAL"), "hid": "0,3"})  # extremo oculto


def test_persistencia_guardar_cerrar_recuperar_repetir(tmp_path):
    path = str(tmp_path / "p.db"); database.init_db(path)
    mod = MODULES[2]; params = defaults(mod.META["fields"], "NORMAL")
    r1 = mod.run("NORMAL", params); sim_id = database.save_simulation(r1, path)
    sim = database.load_simulation(sim_id, path)                 # conexion nueva = "app reiniciada"
    assert sim["params"] == json.loads(json.dumps(params)) and sim["scenario"] == "NORMAL"
    r2 = mod.run(sim["scenario"], sim["params"])
    v1 = {m["name"]: m["value"] for m in r1["methods"]}
    assert all(abs(v1[m["name"]] - m["value"]) < 1e-12 for m in r2["methods"])
    guardados = {e["metodo"]: json.loads(e["resumen"]["valor_json"]) for e in sim["ejecuciones"]}
    assert guardados["Raíz D(t)=D_u · Newton"] == pytest.approx(v1["Raíz D(t)=D_u · Newton"])
    assert any(len(e["iteraciones"]) > 0 for e in sim["ejecuciones"])
    database.delete_simulation(sim_id, path)
    assert database.load_simulation(sim_id, path) is None


def test_rutas_web(client):
    c, _ = client
    assert c.get("/").status_code == 200
    assert c.get("/modulo/3?escenario=MAL").status_code == 200
    assert c.get("/modulo/9").status_code == 404
    r = c.post("/modulo/1", data={"escenario": "NORMAL", "accion": "ejecutar"})
    assert r.status_code == 200 and b"Interpretaci" in r.data
    r = c.post("/modulo/1", data={"escenario": "NORMAL", "accion": "guardar"}, follow_redirects=True)
    assert b"guardada" in r.data
    assert c.get("/historial").status_code == 200 and c.get("/historial/1").status_code == 200
    assert c.get("/historial/777").status_code == 404
    r = c.post("/modulo/1", data={"escenario": "NORMAL", "v0": "-1"})
    assert b"debe ser" in r.data
