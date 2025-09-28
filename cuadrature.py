#!/usr/bin/env python3
r"""
Módulo de cuadratura Gaussiana para integración numérica, usando los polinomios de Legendre.

Este módulo implementa el método de cuadratura Gauss-Legendre para resolver
la integral: I = ∫₁³ [x⁶ - x² sin(2x)] dx


Para esto, se implemento el método mediante un script de python `cuadrature.py`, se probo varios valores de \(N\) con el fin de encontrar un valor para la integral con un error relativo menor a 1e-10, respecto al valor analitico.

"""

import numpy as np
import matplotlib.pyplot as plt


def gaussxw(N):
    r"""
    Calcula los nodos y pesos de la cuadratura Gauss-Legendre en el intervalo [-1, 1].
    
    Utiliza los polinomios de Legendre para generar los puntos de colocación
    y sus respectivos pesos que se utilizaran para la cuadratura Gaussiana.
    
    Args:
        N (int): Número de puntos de la cuadratura, cambiar este valor cambia el resultado de la aproximación de la integral.
        
    Returns:
        pto (numpy.ndarray): Array con los puntos de colocación/nodos
        peso (numpy.ndarray): Array con los pesos correspondientes
            
    Examples:
        >>> pto, peso = gaussxw(1)  # Un solo punto
        >>> pto[0]  # El punto está en el centro
        0.0
        >>> peso[0]  # El peso es 2 (longitud del intervalo [-1,1])
        2.0
        >>> pto, peso = gaussxw(2) # Con 2 puntos
        >>> len(pto)
        2
        >>> len(peso)
        2
    """
    pto, peso = np.polynomial.legendre.leggauss(N)
    return pto, peso


def gaussxwab(liminf, limsup, pto, peso):
    r"""
    Escala los nodos y pesos desde el intervalo [-1, 1] a [liminf, limsup].
    
    Transforma los puntos de colocación y pesos de la cuadratura Gaussiana
    del intervalo estándar [-1, 1] a nuevos puntos de colocación y pesos para cualquier intervalo [a, b].
    
    Args:
        liminf (float): Límite inferior del intervalo de integración.
        limsup (float): Límite superior del intervalo de integración.
        pto (numpy.ndarray): Puntos de colocación en [-1, 1].
        peso (numpy.ndarray): Pesos correspondientes en [-1, 1].
        
    Returns:
        pto_esc (numpy.ndarray): Puntos transformados al nuevo intervalo
        peso_esc (numpy.ndarray): Pesos transformados al nuevo intervalo
            
    Examples:
        >>> pto, peso = gaussxw(2)
        >>> pto_esc, peso_esc = gaussxwab(0, 1, pto, peso)
        >>> len(pto_esc)
        2
        >>> # Los puntos están ahora en [1, 3] en lugar de [-1, 1]
        >>> min(pto_esc) >= 1
        True
        >>> max(pto_esc) <= 3
        True
    """
    return (0.5 * (limsup - liminf) * pto + 0.5 * (limsup + liminf), 
            0.5 * (limsup - liminf) * peso)


def func_arg_int(x):
    r"""
    Función que queremos integrar: f(x) = x⁶ - x² sin(2x).
    
    Args:
        x (float or numpy.ndarray): Variable independiente, preimagen.
        
    Returns:
        (float or numpy.ndarray): Valor de la función evaluada en x.
        
    Examples:
        >>> #Ejemplo con punto
        >>> func_arg_int(0)
        0.0
        >>> func_arg_int(1) > 0
        True
        >>> # Ejemplo con array
        >>> x_vals = np.array([0, 1])
        >>> results = func_arg_int(x_vals)
        >>> len(results) == 2
        True
    """
    return x**6 - x**2 * np.sin(2 * x)


def eva_int(pto, peso, func):
    r"""
    Evalúa la integral aproximada usando cuadratura Gauss-Legendre.
    
    Calcula la integral usando la suma de la función evaluada en los puntos de colocación por los pesos respectivos.
    
    Args:
        pto (numpy.ndarray): Puntos de colocación.
        peso (numpy.ndarray): Pesos respectivos.
        func (callable): Función a integrar.
        
    Returns:
        (float): Valor aproximado de la integral.
        
    Examples:
        >>> pto, peso = gaussxw(2)
        >>> pto_esc, peso_esc = gaussxwab(0, 1, pto, peso)
        >>> # Integrar función constante f(x) = 2
        >>> result = eva_int(pto_esc, peso_esc, lambda x: 2)
        >>> result
        2.0
        >>> # El resultado es de la integral de 2 de 0 a 1 es 2
    """
    return np.sum(peso * func(pto))


def derivada_analitica(x):
    r"""
    Calcula la antiderivada de f(x) = x⁶ - x² sin(2x).
    
    La antiderivada (se calculo manualmente) es: F(x) = x⁷/7 + x²/2 cos(2x) - x/2 sin(2x) - 1/4 cos(2x)
    
    Args:
        x (float or numpy.ndarray): Variable independiente, preimagen.
        
    Returns:
        (float or numpy.ndarray): Valor de la antiderivada evaluada en x.
        
    Examples:
        >>> # F(0) = -1/4
        >>> abs(derivada_analitica(0) - (-0.25)) < 1e-10
        True
        >>> # Verificar que F'(x) = f(x) aproximadamente, usando la definición de derivada y el calculo de un error con tolerancia de 1e-10
        >>> x = 1.0
        >>> h = 1e-10
        >>> derivada_num = (derivada_analitica(x + h) - derivada_analitica(x)) / h
        >>> derivada_exacta = func_arg_int(x)
        >>> abs(derivada_num - derivada_exacta) < 1e-10
        True
    """
    return (x**7 / 7 + x**2 / 2 * np.cos(2*x) - 
            x / 2 * np.sin(2*x) - 1/4 * np.cos(2*x))


def main():
    r"""
    Resuelve la integral definida $\int_{1}^{3} [x^{6} - x^{2} \sin(2x)] dx$ utilizando la Cuadratura Gauss-Legendre y determina el número de puntos (N) necesario para alcanzar una tolerancia de error relativo.

    El script calcula iterativamente la integral para valores crecientes de N
    (número de puntos de cuadratura) hasta que el error relativo de la
    aproximación respecto al valor analítico es menor que una tolerancia
    establecida de $10^{-10}$. Llegamos que con N=7, el error esta debajo de la tolerancia. 

    Outputs:
        - Imprime los resultados de la convergencia en la consola.
        - Genera y guarda dos gráficos PNG:
          1. 'convergencia.png': Gráfico de convergencia de la integral.
          2. 'error.png': Gráfico del error relativo vs N.
    """
    print("Cuadratura Gaussiana: Resolución de ∫₁³ [x⁶ - x² sin(2x)] dx")
    # Calculo del valor analítico/exacto
    analytic = derivada_analitica(3) - derivada_analitica(1)
    print("Valor analítico:", analytic)

    # Se prueba con varios valores de N hasta alcanzar que el error alcanza una tolerencia de 1e-10
    tol = 1e-10
    N = 1
    N_valores = []
    integral_valores = []
    errores = []

    while True:
        # Obtener puntos y pesos de Gauss-Legendre
        pto, peso = gaussxw(N)
        # Escalar al intervalo [1, 3]
        pto_esc, peso_esc = gaussxwab(1, 3, pto, peso)
        # Evaluar integral aproximada
        approx = eva_int(pto_esc, peso_esc, func_arg_int)
        # Calcular error relativo
        err = abs(approx - analytic) / analytic
        
        # Guardar resultados
        N_valores.append(N)
        integral_valores.append(approx)
        errores.append(err)
        
        print(f"N = {N}, integral ≈ {approx:.12f}, error = {err:.12f}")
        
        # Verificar convergencia
        if err < tol:
            break
        N += 1
        # Evitar bucle infinito
        if N > 50:
            break



    print(f"\nPrimer N con error < {tol} es N = {N}, con integral ≈ {approx:.12f}")
    print("Notamos que este valor de N corresponde a lo que esperamos, en efecto, la teoria nos dice que este método es exacto para polinomios hasta grado 2N-1. Aqui se puede aproximar el sen tal que f(x) se aproxima a un polinomio de grado 6.")
    # Gráfico 1: Convergencia
    # Se grafica la integral para varios valores de N para ver que converge al valor teórico
    plt.figure(figsize=(9, 6))
    plt.plot(N_valores, integral_valores, marker='o', label="Aproximación")
    plt.axhline(analytic, color="red", linestyle="--", label="Valor analítico")
    plt.xlabel("N (número de puntos)")
    plt.ylabel("Valor de la integral")
    plt.title("Convergencia de la cuadratura Gauss-Legendre")
    plt.legend()
    plt.grid(True, ls="--", alpha=0.6)
    plt.savefig("convergencia.png", dpi=300)
    print("Gráfico guardado: convergencia.png")


    # Gráfico 2: Error relativo
    # Se grafica el error para diferentes valores de N para ver que de crece
    plt.figure(figsize=(9, 6))
    plt.semilogy(N_valores, errores, marker='o', label="Error relativo")
    plt.xlabel("N (número de puntos)")
    plt.ylabel("Error relativo")
    plt.title("Disminución del error con Gauss-Legendre")
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.legend()
    plt.savefig("error.png", dpi=300)
    print("Gráfico guardado: error.png")

if __name__ == "__main__":
    """
    Script principal que ejecuta el método numérico y encuentra el N para
    obtener un valor de la integral con un error relativo de bajo de la tolerancia.
    """
    main()

