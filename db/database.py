"""Capa de datos SQLite de SIMTELE-RISK (Semana 1: DER inicial). Guardar/recuperar simulaciones se agrega en la Semana 2."""
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
