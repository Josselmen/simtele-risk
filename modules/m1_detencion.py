"""Modulo 1 - Detencion de emergencia de cabinas.

Modelo: x' = v, v' = -(a0 + c*v);  distancia total = v0*tr + x(t_det).
Metodos: RK2 y RK4 (EDO), trapecio/Simpson/Romberg (distancia = integral de v),
biseccion/Newton/secante (v(t) = 0). Escenario MAL CONDICIONADO: raiz doble.

SEMANA 1: solo se plantea el caso (modelo, metodos asignados y definicion de los tres
escenarios). La implementacion numerica se hace en las semanas 2 y 3.
"""

MODULE_ID = 1
META = dict(
    id=1, slug="detencion", nombre="Detención de emergencia", icono="🚡",
    descripcion="¿Cuánto tiempo y distancia necesita una cabina para detenerse tras un frenado de emergencia?",
    modelo="x' = v,  v' = −(a₀ + c·v),  x(0)=0, v(0)=v₀.  Distancia total = v₀·tᵣ + x(t_det).",
    metodos="Pendiente de implementar (Semana 2-3)",
    metodos_futuros="RK2 · RK4 · Trapecio · Simpson · Romberg · Bisección · Newton · Secante",
    implementado=False,
    escenarios={
        "NORMAL": "Frenado nominal con paso h adecuado.",
        "ESTRES": "Mayor velocidad, freno degradado, reacción lenta y pasos grandes.",
        "MAL": "v(t)=v₀(1−t/t_d)² tiene una raíz doble: la bisección no aplica y Newton pierde velocidad.",
    },
    limites=("Modelo unidimensional con parámetros sintéticos: no describe frenos reales, pendiente, viento, "
             "comodidad de pasajeros ni la dinámica del cable."),
    fields=[
        dict(name="v0", label="Velocidad inicial v₀", unit="m/s", type="float", min=0.1, max=30, defaults={"NORMAL": 5.0, "ESTRES": 6.0, "MAL": 5.0}),
        dict(name="a0", label="Desaceleración del freno a₀", unit="m/s²", type="float", min=0.05, max=10, defaults={"NORMAL": 1.0, "ESTRES": 0.6}, only=["NORMAL", "ESTRES"]),
        dict(name="c", label="Resistencia proporcional c", unit="1/s", type="float", min=0, max=5, defaults={"NORMAL": 0.2, "ESTRES": 0.2}, only=["NORMAL", "ESTRES"]),
        dict(name="tr", label="Tiempo de reacción tᵣ", unit="s", type="float", min=0, max=15, defaults={"NORMAL": 1.5, "ESTRES": 3.0}, only=["NORMAL", "ESTRES"]),
        dict(name="Dav", label="Distancia disponible", unit="m", type="float", min=1, max=500, defaults={"NORMAL": 20.0, "ESTRES": 20.0}, only=["NORMAL", "ESTRES"]),
        dict(name="h", label="Paso de RK (h)", unit="s", type="float", min=0.001, max=2, defaults={"NORMAL": 0.1, "ESTRES": 0.5}, only=["NORMAL", "ESTRES"]),
        dict(name="td", label="Tiempo de detención ideal t_d", unit="s", type="float", min=0.5, max=60, defaults={"MAL": 6.0}, only=["MAL"]),
    ],
)
