Se prueban los Types 158 y 534-NoHX
Estanque de altura 2 metros y 4m^3 de volumen (diámetro se deriva de esos datos)
En todos los casos el flujo de agua entrante es 400 L/hr
En todos los casos la temperatura inicial es de 20°C uniforme (aunque es posible setear una temperatura inicial diferente para cada nodo).

Se hacen 4 simulaciones para cada Type:
1: Flujo entra por nodo superior y sale por nodo inferior; temperatura de flujo entrante = 5°C.
2: Flujo entra por nodo superior y sale por nodo inferior; temperatura de flujo entrante = 40°C.
3: Flujo entra por nodo inferior y sale por nodo superior; temperatura de flujo entrante = 5°C.
4: Flujo entra por nodo inferior y sale por nodo superior; temperatura de flujo entrante = 40°C.

Las simulaciones se hacen con "time steps" de 3 minutos y los resultados se registran cada 6 minutos, durante un lapso de 24 horas.
Los resultados se incluyen en carpetas numeradas del 1 al 4 según las diferentes condiciones explicadas más arriba.
Los archivos de resultados se leen de la siguiente forma:
	- Columna 1 corresponde al instante de la toma de datos (en horas)
	- Columna 2 corresponde a la temperatura del flujo que sale del estanque (en °C)
	- Columna 3 corresponde a la temperatura promedio del estanque (en °C)
	- Columnas 4 a 13 corresponden a temperaturas de los nodos 1 a 10 respectivamente  (en °C) (nodo 1 es el superior y nodo 10 es el inferior)


A continuación se detallan algunos parámetros que NO se dejaron en sus valores default:

- En todos los casos el coeficiente de pérdida de calor superior, lateral e inferior es 5 kJ/hr.m^2.K. Este es el valor default para el Type 534 pero no para el 158.
- Las propiedades del agua son las "default" según el Type 534, aunque debieron modificarse en el Type 158:
	- Densidad: 1000 kg/m^3
	- Calor específico: 4.19 kJ/Kg.K
	- Conductividad térmica: 2.14 kJ/hr.m.K
- Por último, en el Type 158 se redujeron el número de "auxiliary heat inputs" y el número de termostatos a cero.

Los demás parámetros e inputs se dejaron en sus valores "default"; estos pueden ser revisados en la carpeta "detalles" incluida junto con cada archivo de resultados.