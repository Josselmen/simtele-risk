PRAGMA foreign_keys = ON;

CREATE TABLE caso_estudio (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    nombre         TEXT NOT NULL UNIQUE,
    modelo_version TEXT NOT NULL
);

CREATE TABLE escenario (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    caso_id        INTEGER NOT NULL REFERENCES caso_estudio(id),
    tipo           TEXT NOT NULL CHECK (tipo IN ('NORMAL','ESTRES','MAL_CONDICIONADO')),
    parametros_json TEXT NOT NULL,
    fuente_datos   TEXT NOT NULL DEFAULT 'sintetico'
);

CREATE TABLE simulacion (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    escenario_id   INTEGER NOT NULL REFERENCES escenario(id),
    fecha          TEXT NOT NULL DEFAULT (datetime('now')),
    usuario        TEXT,
    estado         TEXT NOT NULL DEFAULT 'completada'
);

CREATE TABLE ejecucion_metodo (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    simulacion_id  INTEGER NOT NULL REFERENCES simulacion(id),
    metodo         TEXT NOT NULL,
    tolerancia     REAL,
    paso           REAL,
    max_iter       INTEGER,
    convergio      INTEGER NOT NULL CHECK (convergio IN (0,1)),
    tiempo_ms      REAL
);

CREATE TABLE iteracion_resultado (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    ejecucion_id     INTEGER NOT NULL REFERENCES ejecucion_metodo(id),
    numero           INTEGER NOT NULL,
    aproximacion_json TEXT NOT NULL,
    residual         REAL,
    error_abs        REAL,
    error_rel        REAL
);

CREATE TABLE resultado_resumen (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    ejecucion_id     INTEGER NOT NULL UNIQUE REFERENCES ejecucion_metodo(id),
    referencia       REAL,
    valor_json       TEXT NOT NULL,
    metrica_error_json TEXT,
    interpretacion   TEXT
);

CREATE INDEX idx_escenario_caso ON escenario(caso_id);
CREATE INDEX idx_sim_escenario  ON simulacion(escenario_id);
CREATE INDEX idx_ej_sim         ON ejecucion_metodo(simulacion_id);
CREATE INDEX idx_it_ej          ON iteracion_resultado(ejecucion_id);
