# PPV-LAW

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
forma rigurosa -teniendo en cuenta que los parámetros de la regresión son estimados.
En este caso el resultado se aproxima con una ecuación cuadrática.

![logppv vs logsd ](https://github.com/FGBASTANTE/PPV-LAW/assets/52360383/ac34fb9e-63cf-4bb5-9b8d-3567068805a4)

También determina la carga máxima operante en función de la distancia definido 
un umbral de la vpp, el modelo de la ley de ecala y un nivel de confianza.

![Q vs D](https://github.com/FGBASTANTE/PPV-LAW/assets/52360383/8481d803-7dda-4c76-bb30-a44e79180f27)
