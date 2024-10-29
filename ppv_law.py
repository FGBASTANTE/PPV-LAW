"""
Código para determinar la ley de propagación de vibraciones en el terreno a
partir del registro en voladuras de las velocidades pico de partícula (ppv)
con diferentes cargas y distancias.

Se espera que se haya definido previamente el modelo de distancia escalada
(s_d: Distancia/Carga^beta). Como ejemplo, y típicamente para cargas alargadas:
beta = 1/2, y para cargas esféricas: beta = 1/3.

En el fichero de entrada los valores x son los logaritmos decimales de
las distancias escaladas (log10(s_d)); los valores y son, consecuentemente,
los log10(ppv):
x	y
1.76779	0.2001
0.69139	1.96096
1.55308	1.06786
..............

Se supone modelo de ruido lognormal y ofrece, entre otros resultados, la recta
(curva) de seguridad definido un nivel de confianza nc (nc>0.5 ->50%).

Aquella se calcula tanto de forma aproximada (recta de seguridad) como de
forma rigurosa (intervalo de predicción) -teniendo en cuenta que los parámetros del modelo son estimados.
En este caso el resultado se aproxima con una ecuación cuadrática.

También calcula el intervalo de tolerancia definida la cobertura y el nivel de confianza deseado.
Con el intervalo de tolerancia se "asegura" la cobertura deseada de la población con el nc establecido.

Utilidad con fines docentes

@author: Fernando García Bastante
Universidad de Vigo"""

# se cargan los módulos que se van a emplear
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import scipy.special as sp

''''""""""""""""""""""""""""""""""""""""'''

def ppv_regress(filename, nc=0.9, cobertura=0.95, ngrid=20):
    """
    Calculates the linear regression of the sample (x, y) --> (log_s_d, log_v_e),
    the prediction intervals with a confidence level nc (with nc > 0.5) at the
    points x of the input data, and on a grid of points (x_grid). Plots and saves
    the results, calculates the coverage. The starting data and calculations are
    performed in log10 scale.
    Parameters:
    filename (str): Path to the input file containing 'x' and 'y' columns.
    nc (float): Confidence level for prediction intervals (default is 0.9).
    cobertura (float): Coverage for tolerance intervals (default is 0.95).
    ngrid (int) : Number of points in the grid (default is 20).
    Returns:
    tuple: Contains the following elements:
        - nc_equation_pred (np.ndarray): Coefficients of the quadratic fit for the prediction interval (nc).
        - nc_equation_aprox (tuple): Intercept and slope for the approximate method.
        - nc_equation_tol (np.ndarray): Coefficients of the quadratic fit for the tolerance interval (nc).
    Raises:
    FileNotFoundError: If the input file is not found.
    ValueError: If the input file does not contain 'x' and 'y' columns.
    Notes:
    - The function avoids division by zero by adjusting nc if it is exactly 0.5.
    - The function prints the coverage of the models over the sample.
    - The results are plotted and saved, including the regression line and prediction intervals.
    """

    # para evitar una division entre, prácticamente, 0 si nc=0.5
    if nc == 0.5:
        nc = 0.500001

    # lectura de los datos en fichero de texto con valores (x, y)
    try:
        df_x = pd.read_csv(filename, sep=r'\s+')
        if "x" not in df_x.columns or "y" not in df_x.columns:
            raise ValueError("Input file must contain 'x' and 'y' columns.")
    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")
        return
    except ValueError as e:
        print(f"Error: {e}")
        return

    x = df_x["x"]
    y = df_x["y"]

    # número de datos del registro y el grid de puntos
    n = len(x)
    x_grid = np.linspace(np.min(x), np.max(x), ngrid)

    # regresión lineal, predicción (valor esperado), residuos y su desviación estándar
    slope, intercept = np.polyfit(x, y, deg=1)
    y_predict = intercept + slope * x
    rss_y = y - y_predict
    mse = np.sqrt(np.square(rss_y).sum() / (n - 2))

    # cálculos intermedios para obtener el error total en cada punto
    x_mean = np.mean(x)
    x_gap = x - x_mean
    ss = np.square(x_gap).sum()

    # cálculo del error (en puntos x y x_grid): método aproximado y riguroso
    # (teniendo en cuenta el error del modelo) e intervalo de tolerancia
    se_x_aprox = st.norm.ppf(nc, loc=0, scale=mse)
    se_x = (
        np.sqrt(1 + 1 / n + np.square((x - x_mean)) / ss) *
        mse * st.t.ppf(nc, df=n - 2)
    )
    se_x_grid = (
        np.sqrt(1 + 1 / n + np.square((x_grid - x_mean)) / ss)
        * mse
        * st.t.ppf(nc, df=n - 2)
    )

    # cálculo del intervalo de tolerancia
    _x_tol = np.sqrt(1 / n + np.square((x - x_mean)) / ss)
    _x_tol_grid = np.sqrt(1 / n + np.square((x_grid - x_mean)) / ss)
    zp_d = st.norm.ppf(cobertura) / _x_tol
    zp_d_grid = st.norm.ppf(cobertura) / _x_tol_grid
    se_x_tol = sp.nctdtrit(n - 2, zp_d, nc, out=None) * mse * _x_tol
    se_x_tol_grid = sp.nctdtrit(
        n - 2, zp_d_grid, nc, out=None) * mse * _x_tol_grid

    # cálculo de predicción en x con nivel de confianza nc (ambos métodos e intervalo de tolerancia)
    y_pred_nc_aprox = y_predict + se_x_aprox
    y_pred_nc = y_predict + se_x
    y_tol_nc = y_predict + se_x_tol

    # cálculo de cobertura sobre la muestra con los tres modelos
    cover_aprox = np.mean(y_pred_nc_aprox > y)
    cover_pred = np.mean(y_pred_nc > y)
    cover_tol = np.mean(y_tol_nc > y)
    cover = (cover_aprox, cover_pred, cover_tol)
    print("la cobertura sobre la muestra de los modelos (aprox/pred/toler) es ", cover)

    # se repiten los cálculos para el grid x_grid
    y_predict_x_grid = intercept + slope * x_grid
    y_pred_x_grid_nc_aprox = y_predict_x_grid + se_x_aprox
    y_pred_x_grid_nc = y_predict_x_grid + se_x_grid
    y_tol_x_grid_nc = y_predict_x_grid + se_x_tol_grid

    # ajuste aproximado: recta de seguridad con nc
    nc_equation_aprox = (intercept + se_x_aprox, slope)
    # ajuste del intervalo de predicción a una ecuación de segundo grado
    nc_equation_pred = np.polyfit(x_grid, y_pred_x_grid_nc, deg=2, full="true")
    # ajuste del intervalo de tolerancia a una ecuación de segundo grado
    nc_equation_tol = np.polyfit(x_grid, y_tol_x_grid_nc, deg=2, full="true")

    # gráfico de resultados
    plt.figure(0)
    _ = plt.plot(x, y, marker=".", linestyle="none", label="data")
    plt.plot(x_grid, y_predict_x_grid, linestyle="solid", label="regression")
    plt.plot(x_grid, y_pred_x_grid_nc_aprox,
             linestyle="dashed", label="nc_aprox")
    plt.plot(x_grid, y_pred_x_grid_nc, linestyle="solid", label="prediction")
    plt.plot(x_grid, y_tol_x_grid_nc, linestyle="solid", label="tolerance")
    plt.legend()
    # Label axes
    _ = plt.xlabel("log(sd)")
    _ = plt.ylabel("log(ppv)")
    plt.margins(0.05)

    # guardado de resultados en los puntos x
    df_x["y_pred"] = y_predict
    df_x["y_pred_nc_aprox"] = y_pred_nc_aprox
    df_x["y_pred_nc"] = y_pred_nc
    df_x["y_tol_nc"] = y_tol_nc

    # guardado de resultados en los puntos x_grid
    df_x_grid = pd.DataFrame({"x": x_grid})
    df_x_grid["y_pred"] = y_predict_x_grid
    df_x_grid["y_pred_nc_aprox"] = y_pred_x_grid_nc_aprox
    df_x_grid["y_pred_nc"] = y_pred_x_grid_nc
    df_x_grid["y_tol_nc"] = y_tol_x_grid_nc
    print (df_x_grid)
    return nc_equation_pred[0], nc_equation_aprox, nc_equation_tol[0]


# función auxiliar para crear una tabla de carga operante frente a distancia a
# partir de los resultados de la función: ppv_regress
def cargas_sd(
    nc_equation_pred,
    nc_equation_aprox,
    nc_equation_tol,
    ppvumbral=40,
    beta=0.5,
    d_grid=np.linspace(50, 250, 20),
):
    """
    Calculate the maximum cooperating charges (e.g., kg) as a function of distances
    (d_grid, e.g., m) given a threshold value of the PPV (ppvumbral, e.g., mm/s) and the
    beta value of the scaling law used: s_d: Distance/Load^beta.
    Both the rigorous model and the approximate model are used.
    Parameters:
    nc_equation_pred (list or tuple): Coefficients of the prediction interval quadratic equation.
    nc_equation_aprox (list or tuple): Coefficients of the approximate model's linear equation.
    nc_equation_tol (list or tuple): Coefficients of the tolerance interval quadratic equation.
    Returns:
    pd.DataFrame: DataFrame containing distances (D) and calculated loads (Q).
    """


    # resolución de la ecuación de segundo grado para resolución con intervalo de predicción
    logppv = np.log10(ppvumbral)
    a = nc_equation_pred[0]
    b = nc_equation_pred[1]
    c = nc_equation_pred[2] - logppv
    # raíz menor
    logsd_pred = (-b - np.sqrt((b * b) - 4 * a * c)) / (2 * a)
    sd_pred = 10**logsd_pred
    # cálculo de las cargas para la solución con intervalo de predicción
    cargas_prediccion = np.power(d_grid / sd_pred, 1 / beta)

    # cálculo de las cargas para solución aproximada
    logsd_aprox = (-nc_equation_aprox[0] + logppv) / nc_equation_aprox[1]
    sd_aprox = 10**logsd_aprox
    carga_aprox = np.power(d_grid / sd_aprox, 1 / beta)

    # resolución de la ecuación de segundo grado para resolución con intervalo de tolerancia
    a = nc_equation_tol[0]
    b = nc_equation_tol[1]
    c = nc_equation_tol[2] - logppv
    # raíz menor
    logsd_tol = (-b - np.sqrt((b * b) - 4 * a * c)) / (2 * a)
    sd_tol = 10**logsd_tol
    # cálculo de las cargas para la solución con intervalo de tolerancia
    cargas_tolerancia = np.power(d_grid / sd_tol, 1 / beta)

    # gráfico de resultados
    plt.figure(1)
    _ = plt.plot(d_grid, cargas_prediccion,
                 linestyle="solid", label="Qpred vs D")
    plt.plot(d_grid, carga_aprox, linestyle="dashed", label="Qaprox vs D")
    plt.plot(d_grid, cargas_tolerancia, linestyle="solid", label="Qtol vs D")
    plt.legend()

    # Label axes
    _ = plt.xlabel("Distancia")
    _ = plt.ylabel("Carga")
    plt.margins(0.05)

    # se almacenan los resultados en un dataframe
    df_QvsD = pd.DataFrame({"D": d_grid})
    df_QvsD["Q_prediccion"] = cargas_prediccion
    df_QvsD["Qaprox"] = carga_aprox
    df_QvsD["Qtolerancia"] = cargas_tolerancia

    return df_QvsD


if __name__ == "__main__":

    nc_equation_pred, nc_equation_aprox, nc_equation_tol = ppv_regress(
        "data_ppv_vertical.txt", nc=0.90, cobertura=0.95
    )
    
    tabla_Q_D = cargas_sd(nc_equation_pred, nc_equation_aprox,
                          nc_equation_tol, ppvumbral=50, beta=0.5)
