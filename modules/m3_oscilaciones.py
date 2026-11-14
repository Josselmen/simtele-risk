"""Modulo 3 - Oscilaciones por viento.

Modelo: m x'' + c x' + k x = F0 sin(Om t)  ->  x' = v, v' = (F - c v - k x)/m.
Metodos: RK2 (punto medio) y RK4; trapecio/Simpson/Romberg para I = integral de x^2.
Escenario MAL: sensibilidad a k (problema) frente a crecimiento espurio de RK2 (algoritmo).

SEMANA 1: solo se plantea el caso (modelo, metodos asignados y definicion de los tres
escenarios). La implementacion numerica se hace en las semanas 2 y 3.
"""

MODULE_ID = 3
META = dict(
    id=3, slug="oscilaciones", nombre="Oscilaciones por viento", icono="🌬️",
    descripcion="¿Cómo cambia la oscilación de una cabina ante la fuerza del viento, el amortiguamiento y la rigidez?",
    modelo="m·x'' + c·x' + k·x = F₀·sen(Ω·t);  x' = v,  v' = (F − c·v − k·x)/m.  Exposición I = ∫ x(t)² dt.",
    metodos="Pendiente de implementar (Semana 2-3)",
    metodos_futuros="RK2 · RK4 · Trapecio · Simpson · Romberg",
    implementado=False,
    escenarios={
        "NORMAL": "Vibración libre subamortiguada con solución analítica.",
        "ESTRES": "Viento armónico cerca de la resonancia, menor amortiguamiento, pasos grandes.",
        "MAL": "Casi sin amortiguamiento y T largo: sensibilidad a k y crecimiento espurio de RK2 (dos fenómenos distintos).",
    },
    limites=("Oscilador lineal de un grado de libertad con parámetros sintéticos: no representa el péndulo no lineal de una cabina real, "
             "viento turbulento de banda ancha, interacción con el cable ni amortiguamiento aerodinámico."),
    fields=[
        dict(name="m", label="Masa m", unit="kg", type="float", min=10, max=1e5, defaults={"NORMAL": 2000.0, "ESTRES": 2000.0, "MAL": 2000.0}),
        dict(name="k", label="Rigidez k", unit="N/m", type="float", min=10, max=1e6, defaults={"NORMAL": 8000.0, "ESTRES": 8000.0, "MAL": 8000.0}),
        dict(name="c", label="Amortiguamiento c", unit="N·s/m", type="float", min=0, max=1e5, defaults={"NORMAL": 400.0, "ESTRES": 200.0, "MAL": 16.0}),
        dict(name="F0", label="Amplitud del viento F₀", unit="N", type="float", min=0, max=1e5, defaults={"NORMAL": 0.0, "ESTRES": 600.0, "MAL": 0.0}),
        dict(name="Om", label="Frecuencia del viento Ω", unit="rad/s", type="float", min=0.01, max=50, defaults={"NORMAL": 1.8, "ESTRES": 1.8, "MAL": 1.8}),
        dict(name="x0", label="Posición inicial x₀", unit="m", type="float", min=-10, max=10, defaults={"NORMAL": 0.0, "ESTRES": 0.0, "MAL": 0.0}),
        dict(name="v0", label="Velocidad inicial v₀", unit="m/s", type="float", min=-10, max=10, defaults={"NORMAL": 1.0, "ESTRES": 0.0, "MAL": 1.0}),
        dict(name="T", label="Tiempo total T", unit="s", type="float", min=1, max=1000, defaults={"NORMAL": 20.0, "ESTRES": 60.0, "MAL": 300.0}),
        dict(name="kexp", label="Pasos N = 2^e (exponente e)", unit="", type="int", min=6, max=13, defaults={"NORMAL": 8, "ESTRES": 9, "MAL": 12}),
    ],
)
