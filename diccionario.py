# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 14:05:17 2025

@author: Raul-CDH

"""

#tabla = np.zeros((m, m))  # Matriz vacía con ceros
#dir(np)
#pintar graficas import matplotlib.pyplot as plt

"""VANDERMOD"""
# ejemplo: generar una matriz de tipo Vandermonde  

def vandermonde(nodos): 
    """requiere haber importado numpy as np -
    #nodos es un array con los nodos"""
    n = len(nodos) 
    X = np.zeros((n, n)) 
    for j in range(n): #iteración columnas
        for i in range(n): #iteración filas
            X[i, j] = nodos[j] ** i #coge el punto xj y lo eleva al indice fila
    M = np.asmatrix(X) #convierte a matrix
    return(M) #devuelve matrix

vandermonde([2,3,4])
#####################FACTORIAL#######################################
import math

numero = 5
resultado = math.factorial(numero)
print(resultado)
#DESCOMPONER EN FACTORES PRIMOS
from sympy import factorint

numero = 60
factores = factorint(numero)
print(factores)


#######################################RAICES########################################

"""
    Escribe un procedimiento para calcular las raíces reales de una ecuación de segundo grado que
    permita utilizar las fórmulas adecuadas, según el valor de b, para evitar diferencias entre valores
    próximos, tal y como se ha visto en clase de teoría.
    
"""
import math as mt

def funcion_doblegrado(a, b, c):#ecuacion segundo grado
    """
    hay que meter 3 numeros a,b,c.
    a= será x^2,b = x, c = termino independiente
    creamos una lista vacia para guardar los resultados.
    comprobamos cada cosa, si la a = 0 y b =0 no es una ecuacion
    si a == 0 es eciacion de primer grado.
    
    else: empezamos a clcular el discriminanto (lo que está dentro de la raiz)
    hacemos 3 comprobaciones mas, si discriminante es < 0 no se puede ya que no hay raices negativas
    si > 0 la reiz tiene dos soluciones.
    
    """
    lista_doblegrado = []
    if (a == 0 and b == 0):
        print("No es una ecuacion")
        print(lista_doblegrado)
    elif (a == 0):
        print("Es una ecuacion de primer grado")
    else:
        discriminante = b**2 - 4*a*c
        multi = 2*a
        if(discriminante < 0):
            print("La ecuacion no tiene raices reales")
        elif(discriminante == 0):
            print("Posee raiz doble")
            x=-b/2*a
            lista_doblegrado.append(x)
            print(lista_doblegrado)
        else:
            if(b > 0):
                x1 = multi / (-b - mt.sqrt(discriminante))
                x2 = (-b - mt.sqrt(discriminante)) / multi
                lista_doblegrado.append(x1)
                lista_doblegrado.append(x2)
                print(lista_doblegrado)
            elif(b < 0):
                x1 = (-b + mt.sqrt(discriminante)) / multi
                x2 = multi / (-b + mt.sqrt(discriminante))
                lista_doblegrado.append(x1)
                lista_doblegrado.append(x2)
                print(lista_doblegrado)
                
#(1) x2 + 12345678987x + 1 = 0.
funcion_doblegrado(1, 12345678987, 1)

#(2) x2 − 12345678987x + 1 = 0.
funcion_doblegrado(1, -12345678987, 1)

#(3) x2 − 4x + 4 = 0.
funcion_doblegrado(1, -4, 4)

#(4) x2 − 4x + 5 = 0.
funcion_doblegrado(1, -4, 5)
##############################TYLOR#########################################################
"tenemos que importar todo lo de sympy para poder usar funciones como sin() o symbols("")."
"al poner .remove() estamos eliminando de orden superior dejando solo la expansion"
"parametros: funcion, punto, y orden"
"si se pide aproximar el valor en 3 , remplazamos las x con 3"
from sympy import *
x = symbols("x")
## serie: el polinomio  o serie de Tylor
def serieDeTylor(f, x0, n):
    s = f.series(x, x0, n+1).removeO() 
    return s

serieDeTylor(sin(x),0,5)
###############################LAGRANGE###########################################################
"recordemos que en lagrange lo que buscamos es la forma"
"(x-x1)/(x0-x1) * (x-x2)/(x0-x2)"
"el primer bucle itera por cada elemento o nodo mientras que el bucle interno construya cada base L0"
import numpy as np
x = symbols("x")
nodos = [2, 5/2, 4]
valores_f = [1/nodo for nodo in nodos] 

def polinomio_lagrange(nodos, valores_f):
    """
    se espera recibir una lista de nodos, y valores de la funcion, esdecir (Y)
    """
    grado = len(nodos)
    polinomio = 0  # Inicializamos el polinomio

    for i in range(grado):
        termino_Li = 1  # Base de Lagrange L_i(x)
        for j in range(grado):
            if i != j:
                termino_Li *= ((x - nodos[j]) / (nodos[i] - nodos[j]))
        
        # Sumar término f(x_i) * L_i(x)
        polinomio += valores_f[i] * termino_Li

    return simplify(polinomio)

# Construcción del polinomio de interpolación de Lagrange
polinomio_interpolador = polinomio_lagrange(nodos, valores_f)
print("Polinomio interpolador de Lagrange:")
print(polinomio_interpolador)
help(polinomio_lagrange)
# Evaluación del polinomio en un punto específico
punto = 3  # Por ejemplo, evaluamos en x = 3
evaluacion = polinomio_interpolador.subs(x, punto)

"""otra forma de usar lagrande es con la libreria scipy"""
from scipy import interpolate
import matplotlib.pyplot as plt

# Datos
nodos = [2, 5/2, 4]
valores_f = [1/nodo for nodo in nodos]

# Crear una función interpolante (polinomio de Lagrange implícito en scipy)
f_interpolada = interpolate.lagrange(nodos, valores_f)

# Evaluar la función interpolada en un punto
punto_a_aproximar = 1.5
valor_aproximado = f_interpolada(punto_a_aproximar)

print(f"Valor aproximado de f(1.5) con interpolación de Lagrange: {valor_aproximado}")

###############################NEWTON#####################################################
def polinomi_Interpolador_de_Nweton(nodos,valores_f,x):
    diferencias_divididas = [valores_f.copy()]
    grado_n = len(nodos)
    polinomioDeNewton = 0
    for i in range(1,grado_n): # i, es el orden de la diferencia dividida
        la_fila_actual = [] 
        for j in range(grado_n-i):#recorrer los nodos dsiponibles para el orden j
          numerador = (diferencias_divididas[i-1][j+1] - diferencias_divididas[i-1][j])
          denominador = (nodos[j+i]-nodos[j])
          diferencia = numerador/denominador  
          la_fila_actual.append(diferencia)
        diferencias_divididas.append(la_fila_actual)
        
    polinomioDeNewton = diferencias_divididas[0][0]# El primer término es f(x0)
    factor = 1  # Factor de (x - x0), (x - x1), ...
    
    # Añadimos los términos correspondientes
    for j in range(1, grado_n):
        factor *= (x - nodos[j-1])  # (x - x0), (x - x1), ...
        polinomioDeNewton += diferencias_divididas[j][0] * factor  # Añadimos el término correspondiente al polinomio
    
    return polinomioDeNewton

nodos = [0,1,2]
valores_f = [-2,1,4]
resultado = polinomi_Interpolador_de_Nweton(nodos, valores_f, x)
print(f"Valor aproximado de f({x}) = {resultado}")

# Punto donde queremos evaluar
x = 1.5

# Resultado
resultado = polinomioNweton(nodos, valores_f, x)
print(f"Valor aproximado de f({x}) = {resultado}")

##################################HERMITE################################################
"""
para hallar la derivada de una funcion se usa la funcion de sympy sp.diff(f, x)
"""
import sympy as sp

def polinomioHermite(nodos, valores_f, derivadas_f, x):
    n = len(nodos)
    
    # Construimos una lista de nodos duplicados
    nodos_expandidos = []
    valores_expandidos = []
    
    for i in range(n):
        nodos_expandidos.append(nodos[i])
        valores_expandidos.append(valores_f[i])
        
        nodos_expandidos.append(nodos[i])  # Duplicamos el nodo
        valores_expandidos.append(valores_f[i])  # Misma función
    
    # Matriz de diferencias divididas
    diferencias_divididas = [[0] * (2 * n) for _ in range(2 * n)]
    
    # Llenamos la primera columna con los valores de la función
    for i in range(2 * n):
        diferencias_divididas[i][0] = valores_expandidos[i]
    
    # Segunda columna: si es nodo repetido, usamos derivada
    for i in range(1, 2 * n):
        if nodos_expandidos[i] == nodos_expandidos[i - 1]:  
            diferencias_divididas[i][1] = derivadas_f[i // 2]  # Usamos la derivada
        else:
            diferencias_divididas[i][1] = (diferencias_divididas[i][0] - diferencias_divididas[i-1][0]) / (nodos_expandidos[i] - nodos_expandidos[i-1])
    
    # Llenamos el resto de la tabla de diferencias divididas
    for j in range(2, 2 * n):
        for i in range(2 * n - j):
            diferencias_divididas[i][j] = (diferencias_divididas[i+1][j-1] - diferencias_divididas[i][j-1]) / (nodos_expandidos[i+j] - nodos_expandidos[i])
    
    # Construcción del polinomio
    polinomio = diferencias_divididas[0][0]
    factor = 1
    
    for j in range(1, 2 * n):
        factor *= (x - nodos_expandidos[j-1])
        polinomio += diferencias_divididas[0][j] * factor

    return polinomio

# Ejemplo de uso:
nodos = [0,1]
valores_f = [-2,1]
derivadas_f = [0,0]  # Ejemplo de derivadas
x_aproximar = 1.5

resultado = polinomioHermite(nodos, valores_f, derivadas_f, x_aproximar)
print(f"Valor aproximado de f({x_aproximar}) = {resultado}")

##-------------------------------------------
import numpy as np
import math
def tabla_diferencias_divididas_hermite(nodos, valores_f, derivadas_f):
    n = len(nodos)
    m = 2 * n  # Duplicamos cada nodo en la tabla
    tabla = np.zeros((m, m))  # Matriz vacía con ceros

    # Paso 1: Llenamos la primera columna con los nodos repetidos
    x_expandido = []  # Lista con nodos duplicados
    for i in range(n):
        x_expandido.append(nodos[i])
        x_expandido.append(nodos[i])
    
    # Paso 2: Llenamos la segunda columna con f(x)
    for i in range(n):
        tabla[2 * i][0] = valores_f[i]
        tabla[2 * i + 1][0] = valores_f[i]

    # Paso 3: Llenamos la tercera columna con f'(x) en posiciones pares
    for i in range(n):
        tabla[2 * i + 1][1] = derivadas_f[i]  # Derivada en la segunda columna (posiciones impares)
        if i > 0:
            # Diferencia dividida en la tercera columna (solo posiciones impares)
            tabla[2 * i][1] = (tabla[2 * i][0] - tabla[2 * i - 1][0]) / (x_expandido[2 * i] - x_expandido[2 * i - 1])

    # Paso 4: Llenamos el resto de la tabla
    for j in range(2, m):  # Columnas de la tabla
        for i in range(m - j):  # Filas
            tabla[i][j] = (tabla[i + 1][j - 1] - tabla[i][j - 1]) / (x_expandido[i + j] - x_expandido[i])

    return x_expandido, tabla

# Ejemplo de uso:
nodos = [0, 1]
valores_f = [-2, 1]
derivadas_f = [0, 0]  # Derivadas en cada nodo

x_expandido, tabla = tabla_diferencias_divididas_hermite(nodos, valores_f, derivadas_f)

# Imprimir la tabla
import pandas as pd
df = pd.DataFrame(tabla, index=x_expandido)
print(df)

#############################COEFINIENTES INDETERMINADOS Integracion########################

import numpy as np
import sympy as sp
import math
def tres_octavos(f, a, b):
    # Definimos los nodos
    x0 = a
    x1 = (2*a + b) / 3
    x2 = (a + 2*b) / 3
    x3 = b
    
    # Planteamos el sistema de ecuaciones para los coeficientes
    A = np.array([
        [1, 1, 1, 1],
        [x0, x1, x2, x3],
        [x0**2, x1**2, x2**2, x3**2],
        [x0**3, x1**3, x2**3, x3**3]
    ])
    
    B = np.array([
        b - a,
        (b**2 - a**2) / 2,
        (b**3 - a**3) / 3,
        (b**4 - a**4) / 4
    ])
    
    # Resolvemos el sistema de ecuaciones
    coeficientes = np.linalg.solve(A, B)
    
    # Definimos la fórmula de tres octavos
    integral = (coeficientes[0] * f(x0) + coeficientes[1] * f(x1) + coeficientes[2] * f(x2) + coeficientes[3] * f(x3))
    
    return integral
# Ejemplo de uso
def f(x):
    return x

def ff(x):
    return sp.cos(x)

def fff(x):
    return sp.exp(x)
resultado = tres_octavos(f, 0, 1)
resultado2 = tres_octavos(ff, -math.pi/2, math.pi/2)
resultado3 = tres_octavos(fff, 0, 1)
print(f"La aproximación de la integral definida de x en el intervalo [0, 1] es: {resultado:.15f}")
print(f"La aproximación de la integral definida de cos(x) en el intervalo [-pi/2, pi/2] es: {resultado2:.15f}")
print(f"La aproximación de la integral definida de exp(x) en el intervalo [0, 1] es: {resultado3:.15f}")


###############################COEFICIENTES INDETERMINADOS DERIVACION#############
import numpy as np

def derivada_fpc(f, n, x0, h):
    if n == 1:
        # Primera derivada cinco puntos
        derivada = (-f(x0 + 2*h) + 8*f(x0 + h) - 8*f(x0 - h) + f(x0 - 2*h)) / (12*h)
    elif n == 2:
        # Segunda derivada cinco puntos
        derivada = (-f(x0 + 2*h) + 16*f(x0 + h) - 30*f(x0) + 16*f(x0 - h) - f(x0 - 2*h)) / (12*(h**2))
    elif n == 3:
        # Tercera derivada cinco puntos
        derivada = (f(x0 + 2*h) - 2*f(x0 + h) + 2*f(x0 - h) - f(x0 - 2*h)) / (2*(h**3))
    return derivada

# f(x) = sin(x)
def f(x):
    return np.sin(x)
x0 = 0 # Punto de evaluaciÃ³n
h = 0.1 # Incremento

primera_derivada = derivada_fpc(f, 1, x0, h)
segunda_derivada = derivada_fpc(f, 2, x0, h)
tercera_derivada = derivada_fpc(f, 3, x0, h)

print(f"Primera derivada (n=1): {primera_derivada:.15f}")
print(f"Segunda derivada (n=2): {segunda_derivada:.30f}")
print(f"Tercera derivada (n=3): {tercera_derivada:.15f}")


################################SIMPSION################################################
"""
EJERCICIO 4

Escribe procedimientos en python para implementar las fórmulas de cuadratura
compuestas (rectángulo, punto medio, trapecio y Simpson) con n+1 nodos, y aplíquense para
aproximar el valor de la integral

INTEGRAL [2,0] X^2e^-X^2 dx
"""
import numpy as np
import scipy.integrate as spi
import sympy as sp

# Definición de la función
def f(x):
    return x**2 * np.exp(-x**2)

# Procedimiento para la fórmula del rectángulo compuesto
def rectangulo_compuesto(f, a, b, n):
    h = (b - a) / n
    integral = 0
    for i in range(n):
        integral += f(a + i*h) * h
    return integral

# Procedimiento para la fórmula del punto medio compuesto
def punto_medio_compuesto(f, a, b, n):
    h = (b - a) / n
    integral = 0
    for i in range(n):
        integral += f(a + (i + 0.5)*h) * h
    return integral

# Procedimiento para la fórmula del trapecio compuesto
def trapecio_compuesto(f, a, b, n):
    h = (b - a) / n
    integral = (f(a) + f(b)) / 2
    for i in range(1, n):
        integral += f(a + i*h)
    integral *= h
    return integral

# Procedimiento para la fórmula de Simpson compuesta
def simpson_compuesto(f, a, b, n):
    if n % 2 == 1:
        raise ValueError("El número de subintervalos debe ser par para la fórmula de Simpson.")
    h = (b - a) / n
    integral = f(a) + f(b)
    for i in range(1, n, 2):
        integral += 4 * f(a + i*h)
    for i in range(2, n-1, 2):
        integral += 2 * f(a + i*h)
    integral *= h / 3
    return integral

# Parámetros del problema
a = 0
b = 2
n = 8

# Cálculo de las aproximaciones con las fórmulas compuestas
rectangulo_resultado = rectangulo_compuesto(f, a, b, n)
punto_medio_resultado = punto_medio_compuesto(f, a, b, n)
trapecio_resultado = trapecio_compuesto(f, a, b, n)
simpson_resultado = simpson_compuesto(f, a, b, n)

# Cálculo del valor exacto con sympy
x = sp.Symbol('x')
integral_exacta = sp.integrate(x**2 * sp.exp(-x**2), (x, a, b))

# Cálculo del valor aproximado con scipy
integral_scipy, _ = spi.quad(f, a, b)

# Resultados
print(f"Aproximación con la fórmula del rectángulo compuesto: {rectangulo_resultado:.15f}")
print(f"Aproximación con la fórmula del punto medio compuesto: {punto_medio_resultado:.15f}")
print(f"Aproximación con la fórmula del trapecio compuesto: {trapecio_resultado:.15f}")
print(f"Aproximación con la fórmula de Simpson compuesta: {simpson_resultado:.15f}")
print(f"Valor exacto calculado con sympy: {integral_exacta.evalf()}")
print(f"Valor aproximado calculado con scipy: {integral_scipy}")

#################################################################################
import numpy as np
import math
from math import sin, factorial

def derivacion_cincoPuntos(f, n, x0, h):
    """
    Aproxima la derivada de orden n (n=1,2,3) de la función f en x0
    usando la fórmula central de 5 puntos con coeficientes indeterminados.
    
    la formula de la primera derivada:
        siendo h la distancia entre puntos
        f(x) = −f(xi−2)+8f(xi−1)−8f(xi+1)+f(xi+2)/12h
​

    """
    if n not in [1, 2, 3]:
        raise ValueError("El orden de derivación n debe ser 1, 2 o 3.")

    # Desplazamientos
    d = np.array([-2, -1, 0, 1, 2], dtype=float)
    nodes = x0 + h*d

    # Matriz de Vandermonde
    V = np.array([[d_i**j for d_i in d] for j in range(5)])

    # Vector de la derecha
    b = np.zeros(5)
    b[n] = math.factorial(n) / (h**n)

    # Coeficientes A
    A = np.linalg.solve(V, b)

    # f(nodes) con list comprehension si f es escalar
    fvals = [f(val) for val in nodes]

    # Producto punto
    deriv_approx = np.dot(A, fvals)
    return deriv_approx

# Probamos
x0 = 0
h = 0.1

deriv1 = derivacion_cincoPuntos(sin, 1, x0, h)
deriv2 = derivacion_cincoPuntos(sin, 2, x0, h)
deriv3 = derivacion_cincoPuntos(sin, 3, x0, h)

print("Primera derivada (n=1) =", deriv1)
print(f"Segunda derivada (n=2) = {deriv2:.15f}")
print("Tercera derivada (n=3) =", deriv3)


