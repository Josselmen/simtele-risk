"""Modulos de casos de estudio de SIMTELE-RISK."""
from . import m1_detencion, m2_deterioro, m3_oscilaciones, m4_sensores

MODULES = {1: m1_detencion, 2: m2_deterioro, 3: m3_oscilaciones, 4: m4_sensores}

# Relacion entre modelos y familias de metodos (planteamiento de la Semana 1).
RELACION_HEADERS = ["Familia de métodos", "M1 Detención", "M2 Deterioro", "M3 Oscilaciones", "M4 Sensores", "Semana"]
RELACION = [
    ["Errores (transversal)", "Sí", "Sí", "Sí", "Sí", "1"],
    ["Sistemas lineales (Jacobi, Gauss-Seidel, LU, CG)", "—", "Calibración SPD", "Auxiliar (opcional)", "LU en spline", "1"],
    ["Raíces (bisección, Newton, secante)", "v(t) = 0", "D(t) − D_u = 0", "—", "Cruce de umbral", "2"],
    ["Interpolación (Lagrange, Newton, spline)", "—", "—", "—", "Sí", "2"],
    ["Integración (trapecio, Simpson, Romberg)", "x = ∫v dt", "Tasa (opcional)", "I = ∫x² dt", "—", "3"],
    ["EDO (RK2, RK4)", "(x, v)", "—", "(x, v)", "—", "3"],
]
