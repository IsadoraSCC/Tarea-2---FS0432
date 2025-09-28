# Introducción

El objetivo de esta tarea es resolver la siguiente integral:

\[
I = \int_1^3 \left( x^6 - x^2 \sin(2x) \right)\, dx
\]

Para esto, utilizamos el método de **la cuadratura Gaussiana**, usando los **polinomios de Legendre**, método tambien conocido como **cuadratura Gauss-Legendre**.

Para esto, se implemento el método mediante un script de python `cuadrature.py`, se probo varios valores de \(N\) con el fin de encontrar un valor para la integral con un error relativo menor a $10^(10)$, respecto al valor analitico.

Ademas, se realiza la documentacion del codigo usando **MkDocs**  

