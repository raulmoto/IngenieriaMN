# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 14:05:17 2025

@author: Raul-CDH

"""

#tabla = np.zeros((m, m))  # Matriz vacía con ceros

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

