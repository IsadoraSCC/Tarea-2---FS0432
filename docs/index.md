# Introducción



El objetivo de esta tarea es resolver la siguiente integral:

$$
I = \int_1^3 \left( x^6 - x^2 \sin(2x) \right)\, dx
$$

Para esto, utilizamos el método de **la cuadratura Gaussiana**, usando los **polinomios de Legendre**, método también conocido como **cuadratura Gauss-Legendre**.

Se implementó el método mediante un script de Python `cuadrature.py`.  
Se probaron varios valores de $N$ con el fin de encontrar un valor para la integral con un error relativo menor a $10^{-10}$, respecto al valor analítico.

Además, se realiza la documentación del código usando **MkDocs**.

