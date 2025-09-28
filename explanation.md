# Explicación del Método de Cuadratura Gaussiana

La **cuadratura Gauss-Legendre** es un método que permite aproximar integrales en el intervalo \([-1,1]\):

\[
\int_{-1}^1 f(t)\,dt \approx \sum_{i=1}^N w_i f(t_i)
\]

donde \(t_i\) son las raíces del polinomio de Legendre de orden \(N\), y \(w_i\) son pesos asociados.

Este método permite obtener un mejor resultado que el método de trapecio para un mismo \(N\). 

Se utiliza los polinomios de Legrendre, los cuales son funciones matemáticas que se usan como "base" para nuestro método. Los primeros son:

- \(P_0(x) = 1\)
- \(P_1(x) = x\)
- \(P_2(x) = \frac{1}{2}(3x^2 - 1)\)
- \(P_3(x) = \frac{1}{2}(5x^3 - 3x)\)


Para integrar en un intervalo \([a,b]\), se usa el cambio de variable:

\[
x = \tfrac{b-a}{2} t + \tfrac{b+a}{2}
\]

Esto nos da: 

\[
\int_a^b f(x)\,dx \approx \tfrac{b-a}{2} \sum_{i=1}^N w_i \, f\!\left(\tfrac{b-a}{2}t_i + \tfrac{b+a}{2}\right).
\]


Para nuestro caso específico \([1, 3]\): 
$$x = \frac{3-1}{2}\xi + \frac{1+3}{2} = \xi + 2$$

El algoritmo para realizar la aproximación sigue los siguientes pasos : 
1. Empezar con N = 1 punto
2. Calcular puntos y pesos en [-1, 1] usando Legendre
3. Escalar al intervalo [1, 3]
4. Evaluar la integral usando la suma ponderada
5. Calcular el error comparando con el valor exacto
6. Si el error es muy grande, aumentar N y repetir


Este método es exacto para polinomios hasta grado \(2N-1\).

