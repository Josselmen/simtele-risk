"""Modulo 4 - Reconstruccion de senales de sensores.

Senal sintetica y = sin t en [0, 2*pi]; se ocultan nodos (registros faltantes) y se
reconstruyen con Lagrange, Newton (diferencias divididas) y spline cubico natural.
Raices: cruce del umbral y = yu sobre el spline.

SEMANA 1: solo se plantea el caso (modelo, metodos asignados y definicion de los tres
escenarios). La implementacion numerica se hace en las semanas 2 y 3.
"""

MODULE_ID = 4
META = dict(
    id=4, slug="sensores", nombre="Señales de sensores", icono="📡",
    descripcion="¿Cómo reconstruir registros faltantes de un sensor y cuánto distorsiona el ruido la reconstrucción?",
    modelo="y(t) = sen t en [0, 2π], nodos equiespaciados; yᵢ = sen tᵢ + σ·ξᵢ (ξ ~ N(0,1)); se ocultan algunos nodos.",
    metodos="Pendiente de implementar (Semana 2-3)",
    metodos_futuros="Lagrange · Newton (dif. divididas) · Spline cúbico · Bisección · Newton · Secante",
    implementado=False,
    escenarios={
        "NORMAL": "11 nodos, 2 ocultos, sin ruido.",
        "ESTRES": "15 nodos, 4 ocultos y ruido σ = 0.05; se estudia el efecto del nivel de ruido.",
        "MAL": "25 nodos equiespaciados con ruido 1e-3: la constante de Lebesgue amplifica el ruido en la interpolación polinómica.",
    },
    limites=("Señal sintética y ruido gaussiano independiente: no representa sensores reales con saltos o ruido de otro tipo. "
             "No se extrapola fuera del intervalo muestreado ni se infiere la causa de la falla."),
    fields=[
        dict(name="n", label="Número de nodos n", unit="", type="int", min=7, max=40, defaults={"NORMAL": 11, "ESTRES": 15, "MAL": 25}),
        dict(name="hid", label="Índices ocultos (separados por coma)", unit="", type="text", defaults={"NORMAL": "3,7", "ESTRES": "2,5,8,11", "MAL": "6,12,18"}),
        dict(name="sigma", label="Ruido σ", unit="", type="float", min=0, max=0.5, defaults={"NORMAL": 0.0, "ESTRES": 0.05, "MAL": 0.001}),
        dict(name="seed", label="Semilla aleatoria", unit="", type="int", min=0, max=10 ** 6, defaults={"NORMAL": 42, "ESTRES": 42, "MAL": 11}),
        dict(name="trials", label="Realizaciones para estadísticas", unit="", type="int", min=1, max=200, defaults={"NORMAL": 30, "ESTRES": 50, "MAL": 30}),
        dict(name="yu", label="Umbral de cruce y_u", unit="", type="float", min=0.1, max=0.9, defaults={"NORMAL": 0.5, "ESTRES": 0.5, "MAL": 0.5}),
    ],
)
