# SIMTELE-RISK

Simulador **educativo** de métodos numéricos aplicado a cuatro problemas de transporte por cable.
Compara métodos implementados desde cero en tres escenarios: **NORMAL**, **ESTRÉS** y **MAL CONDICIONADO**.

> Proyecto académico. **Todos los datos y parámetros son sintéticos**; no representan equipos ni datos reales.

## Módulos

| # | Caso de estudio | Métodos numéricos |
|---|-----------------|-------------------|
| 1 | Detención de emergencia | RK2, RK4, trapecio, Simpson, Romberg, bisección, Newton, secante |
| 2 | Deterioro y mantenimiento | Bisección, Newton, secante, Jacobi, Gauss-Seidel, LU, gradiente conjugado |
| 3 | Oscilaciones por viento | RK2, RK4, trapecio, Simpson, Romberg |
| 4 | Señales de sensores | Lagrange, Newton (dif. divididas), spline cúbico, bisección, Newton, secante |

Cada módulo muestra: parámetros validados, tablas de resultados, errores (absoluto, relativo, residuo, orden observado),
tablas de iteraciones, 2 o más gráficos, y **dos conclusiones separadas** (numérica y del problema). Las simulaciones se
pueden **guardar** en SQLite, abrir desde el **Historial**, relanzar y eliminar.

## Requisitos

Python 3.10 o superior.

## Instalación y ejecución

```bash
git clone <URL-DEL-REPOSITORIO>
cd simtele-risk
python -m venv .venv
# Windows:  .venv\Scripts\activate        Linux/Mac:  source .venv/bin/activate
pip install -r requirements.txt
python app.py
```

Abrir <http://127.0.0.1:5000> en el navegador. La base de datos `db/simtele.db` se crea sola al iniciar.

## Pruebas

```bash
python -m pytest -q
```

Verifican cada método contra soluciones analíticas (√2, ∫sen = 2, y' = −y), el orden de convergencia, el comportamiento
en casos mal condicionados, los 12 módulos×escenarios, la validación de parámetros, la persistencia
(guardar → recuperar → repetir) y las rutas web.

## Estructura

```
simtele-risk/
├── app.py                  # aplicación Flask (rutas)
├── requirements.txt
├── README.md
├── db/
│   ├── schema.sql          # 6 tablas: caso_estudio, escenario, simulacion, ejecucion_metodo, iteracion_resultado, resultado_resumen
│   └── database.py         # crear, guardar, listar, recuperar y eliminar simulaciones
├── numerics/               # SOLO matemática (sin Flask)
│   ├── errors.py           # errores, residuo, amplificación
│   ├── linear.py           # LU, Jacobi, Gauss-Seidel, gradiente conjugado
│   ├── interpolation.py    # Lagrange, Newton, spline cúbico natural, constante de Lebesgue
│   ├── roots.py            # bisección, Newton, secante
│   ├── integration.py      # trapecio, Simpson, Romberg
│   └── ode.py              # RK2 (punto medio), RK4
├── modules/                # un archivo por caso de estudio
│   ├── common.py           # validación, tablas, gráficos
│   ├── m1_detencion.py
│   ├── m2_deterioro.py
│   ├── m3_oscilaciones.py
│   └── m4_sensores.py
├── templates/              # páginas HTML (Jinja2)
├── static/style.css
├── tests/                  # pruebas con pytest
└── docs/                   # informes de las semanas 1, 2 y 3
```

## Escenarios

- **NORMAL:** configuración nominal con referencia analítica.
- **ESTRÉS:** parámetros extremos, pasos grandes o ruido.
- **MAL CONDICIONADO:** subexperimento numérico con evidencia cuantitativa (raíz doble, matriz de Hilbert,
  sensibilidad a la rigidez, constante de Lebesgue). Se distingue de un fallo de convergencia.

## Notas técnicas

- Normas: errores y residuos en norma infinito; número de condición en norma 2.
- Iterativos: parada con residuo relativo < tolerancia; Jacobi/Gauss-Seidel hasta 20 000 iteraciones.
- En el Módulo 3 con fuerza externa, la referencia es RK4 con 32 768 pasos (en los informes se usó DOP853 de SciPy).

## Cómo subirlo a GitHub

```bash
git init
git add .
git commit -m "SIMTELE-RISK: version inicial"
git branch -M main
git remote add origin https://github.com/<usuario>/simtele-risk.git
git push -u origin main
```

(Antes, crear en github.com un repositorio vacío llamado `simtele-risk`.) Conviene que **cada integrante haga sus propios commits**.

## Integrantes

[completar]
