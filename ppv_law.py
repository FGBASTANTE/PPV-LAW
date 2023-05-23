# -*- coding: utf-8 -*-
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

Aquella se calcula tanto de de forma aproximada (recta de seguridad) como de
forma rigurosa -teniendo en cuenta que los parámetros del modelo son estimados.
En este caso el resultado se aproxima con una ecuación cuadrática.

@author: Fernando García Bastante
Universidad de Vigo"""

# se cargan los módulos que se van a emplear
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

''''""""""""""""""""""""""""""""""""""""''' 

def ppv_regress(filename, nc=0.9):
    '''Calcula la regresión lineal de la muestra (x, y)--> (log_s_d, log_v_e),
    los intervalos de predicción con un nivel de confianza nc (con nc>0.5) en
    los puntos x de los datos de entrada, y en un grid de puntos (x_grid).
    Grafica y guarda los resultados, calcula la cobertura...
    Los datos de partida y cálculos realizados en escala log10'''
    
    # para evitar una division entre, prácticamente, 0 si nc=0.5
    if nc == 0.5:nc=0.5001
    
    # lectura de los datos en fichero de texto con valores (x, y)
    df_x = pd.read_csv(filename, sep="\s+")    
    x = df_x['x']
    y = df_x['y']
    
    # número de datos de trabajo y el grid de puntos (20 puntos por ejemplo)
    n = len(x)
    x_grid = np.linspace(np.min(x), np.max(x), 20)
    
    # regresión lineal, predicción, residuos y su desviación estándar
    slope, intercept = np.polyfit(x, y, deg=1)
    y_predict = intercept + slope*x
    rss_y = y - y_predict
    mse = np.sqrt(np.square(rss_y).sum()/(n-2))

    # cálculos intermedios para obtener el error total en cada punto
    x_mean = np.mean(x)
    x_gap = x-x_mean
    ss = np.square(x_gap).sum()
    
    # cálculo del error (en puntos x y x_grid): método aproximado y riguroso
    # (teniendo en cuenta el error del modelo)
    se_x_aprox = st.norm.ppf(nc, loc = 0, scale = mse)
    se_x = np.sqrt(1 + 1/n + np.square((x - x_mean))/ss)* \
            mse*st.t.ppf(nc, df = n-2)    
    se_x_grid = np.sqrt(1 + 1/n + np.square((x_grid - x_mean))/ss)* \
            mse*st.t.ppf(nc, df = n-2)
    
    # cáculo de predicción en x con nivel de confianza nc (ambos métodos)
    y_pred_nc_aprox = y_predict + se_x_aprox   
    y_pred_nc = y_predict + se_x    
    
    # cálculo de cobertura con ambos métodos
    cover_aprox = np.mean(y_pred_nc_aprox > y)
    cover = np.mean(y_pred_nc > y)    
    cover = (cover_aprox, cover)
    print('la cobertura de los modelos es ', cover)
    
    # se repiten los cálculos para el grid x_grid
    y_predict_x_grid = intercept + slope*x_grid
    y_pred_x_grid_nc_aprox =  y_predict_x_grid + se_x_aprox
    y_pred_x_grid_nc =  y_predict_x_grid + se_x_grid    
    
    # ajuste aproximado: recta se deguridad con nc
    nc_equation_aprox = (intercept+se_x_aprox, slope)
    # ajuste del método riguroso a una ecuación de segundo grado
    nc_equation = np.polyfit(x_grid, y_pred_x_grid_nc, deg=2, full='true')
    
    # gráfico de resultados
    plt.figure(0)
    _ = plt.plot(x, y, marker='.', linestyle='none', label='data')
    plt.plot(x_grid, y_predict_x_grid, linestyle='solid', label='regression')
    plt.plot(x_grid, y_pred_x_grid_nc_aprox, linestyle='dashed', label='nc_aprox')
    plt.plot(x_grid, y_pred_x_grid_nc, linestyle='solid', label='nc')
    plt.legend()
    # Label axes
    _ = plt.xlabel('log(sd)')
    _ = plt.ylabel('log(ppv)')
    plt.margins(0.05)
    
    # guardado de resultados en los puntos x
    df_x['y_pred'] = y_predict
    df_x['y_pred_nc_aprox'] = y_pred_nc_aprox
    df_x['y_pred_nc'] = y_pred_nc
    
    # guardado de resultados en los puntos x_grid
    df_x_grid = pd.DataFrame({'x': x_grid})
    df_x_grid['y_pred'] = y_predict_x_grid
    df_x_grid['y_pred_nc_aprox'] = y_pred_x_grid_nc_aprox
    df_x_grid['y_pred_nc'] = y_pred_x_grid_nc
    
    return  intercept, slope, df_x, df_x_grid, nc_equation, nc_equation_aprox

results = ppv_regress("data_ppv_vertical.txt", nc=0.95)

# función auxiliar para crear una tabla de carga operante frente a distancia a
# partir de los resultados de la función: ppv_regress
def cargas_sd(nc_equation, nc_equation_aprox):
    ''' cálculo de las cargas máximas operantes (vg. kg) en función de las 
    distancias (d_grid, vg. m) definido un valor umbral de la ppv (vg. mm/s) y
    el valor de beta de la ley de escala empleda: s_d: Distancia/Carga^beta.
    Se emplean el modelo riguroso y el aproximado'''
    
    # valores definidos por el usuario (unidades en ejemplo: mm/s, m y kg)
    ppvumbral=40 # mm/s
    beta=0.5
    d_grid = np.linspace(50, 250, 20)  # m
    
    # resolución de la ecuación de segundo grado para resolución rigurosa
    logppv= np.log10(ppvumbral)
    a = nc_equation[0]
    b = nc_equation[1]
    c = nc_equation[2]-logppv
    
    # raíz menor
    logsd = (-b - np.sqrt((b * b) - 4 * a * c))/(2 * a)
    sd = 10**logsd
    
    # cálculo de las cargas
    cargas = np.power(d_grid/sd, 1/beta)
    
    # cálculo de las cargas para solución aproximada
    logsd_aprox = (-nc_equation_aprox[0]+logppv)/nc_equation_aprox[1]
    sd_aprox = 10**logsd_aprox
    carga_aprox = np.power(d_grid/sd_aprox, 1/beta)
    
    # gráfico de resultados
    plt.figure(1)
    _ = plt.plot(d_grid, cargas, linestyle='solid', label='Q vs D')
    plt.plot(d_grid, carga_aprox, linestyle='dashed', label='Qaprox vs D')
    plt.legend()
    # Label axes
    _ = plt.xlabel('Distancia')
    _ = plt.ylabel('Carga')
    plt.margins(0.05)
    
    # se almacenan los resultados en un dataframe
    df_QvsD = pd.DataFrame({'D': d_grid})
    df_QvsD['Q'] = cargas
    df_QvsD['Qaprox'] = carga_aprox
    return df_QvsD

tabla_Q_D = cargas_sd(results[4][0], results[5])
