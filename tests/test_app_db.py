"""Semana 1: modulo 2 (3 escenarios), modulos planteados, DER y rutas web."""
import sqlite3

import pytest

from app import create_app
from db import database
from modules import MODULES
from modules.common import ValidationError, defaults, parse_params


@pytest.fixture()
def client(tmp_path):
    app = create_app(str(tmp_path / "test.db"))
    app.config["TESTING"] = True
    return app.test_client()


@pytest.mark.parametrize("sc", ["NORMAL", "ESTRES", "MAL"])
def test_modulo2_tres_escenarios(sc):
    mod = MODULES[2]
    r = mod.run(sc, defaults(mod.META["fields"], sc))
    assert r["tables"] and r["charts"] and r["methods"]
    assert r["interpretation"]["numeric"] and r["interpretation"]["problem"]


def test_cuatro_modulos_con_tres_escenarios_definidos():
    assert set(MODULES) == {1, 2, 3, 4}
    for mod in MODULES.values():
        assert set(mod.META["escenarios"]) == {"NORMAL", "ESTRES", "MAL"}
    assert [m.META["implementado"] for m in MODULES.values()] == [False, True, False, False]


def test_validacion_de_parametros():
    _, errs = parse_params(MODULES[2].META["fields"], {"D0": "-3", "n": "abc"}, "NORMAL")
    assert len(errs) == 2
    with pytest.raises(ValidationError):
        MODULES[2].run("NORMAL", {**defaults(MODULES[2].META["fields"], "NORMAL"), "D0": 80.0})   # D0 >= Du


def test_der_tablas_y_casos(tmp_path):
    path = str(tmp_path / "d.db"); database.init_db(path)
    con = sqlite3.connect(path)
    tablas = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    assert {"caso_estudio", "escenario", "simulacion", "ejecucion_metodo", "iteracion_resultado", "resultado_resumen"} <= tablas
    assert con.execute("SELECT COUNT(*) FROM caso_estudio").fetchone()[0] == 4


def test_rutas_web(client):
    c = client
    r = c.get("/")
    assert r.status_code == 200 and b"Relaci" in r.data
    for mid in (1, 2, 3, 4):
        assert c.get(f"/modulo/{mid}?escenario=MAL").status_code == 200
    assert c.get("/modulo/9").status_code == 404
    r = c.post("/modulo/2", data={"escenario": "NORMAL", "accion": "ejecutar"})
    assert r.status_code == 200 and b"Interpretaci" in r.data
    r = c.post("/modulo/2", data={"escenario": "NORMAL", "D0": "-1"})
    assert b"debe ser" in r.data
    assert c.get("/historial").status_code == 404
