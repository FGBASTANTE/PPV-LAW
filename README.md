# PPV-Attenuation_Law

Código para determinar la ley de propagación de vibraciones en el terreno a
partir del registro en voladuras de las velocidades pico de partícula (ppv)
con diferentes cargas y distancias.

Se espera que se haya definido previamente el modelo de distancia escalada
(s_d: Distancia/Carga^beta). Como ejemplo, y típicamente para cargas alargadas:
beta = 1/2, y para cargas esféricas: beta = 1/3.

En el fichero de entrada los valores x son los logaritmos decimales de
las distancias escaladas (log10(s_d)); los valores y son, consecuentemente,
los valores log10(ppv):\
x	y\
1.76779	0.2001\
0.69139	1.96096\
1.55308	1.06786\
..............

Las unidades las define el usuario.

Se supone modelo de ruido lognormal y ofrece, entre otros resultados, la recta
(curva) de seguridad definido un nivel de confianza nc (nc>0.5 ->50%).

Aquella se calcula tanto de forma aproximada (recta de seguridad) como de
forma rigurosa -teniendo en cuenta que los parámetros de la regresión son estimados (intervalo de predicción). En este caso el resultado se aproxima con una ecuación cuadrática.

![image](https://github.com/user-attachments/assets/a8555da2-88e2-479f-aa7d-b313a578ed62)

También calcula el intervalo de tolerancia definida una cobertura de la población y el nc.

Finalmente, determina la carga máxima operante en función de la distancia definido 
un umbral de la ppv, el modelo de la ley de escala y un nivel de confianza, y la cobertura.

![image](https://github.com/user-attachments/assets/6c263f17-3a58-492e-8977-b080497aae72)

