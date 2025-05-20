# -*- coding: utf-8 -*-
"""
    Ejercicio 1 (Ejemplo 1, Tema 4).
    
    Resuelve la ecuación f (x) = x^3 + 4x^2 − 10 = 0 en el intervalo [1, 2]
    utilizando el método de bisección. Para ello, escribe un procedimiento en python que tome como
    parámetros una función f y los extremos a, b de un intervalo.
    Asimismo, representa la función anterior en el intervalo adecuado para comprobar si se verifican las
    hipótesis del teorema de Bolzano.
"""



########Nota
"""
    Si f(c)= 0 encontraste la raíz exacta.
    
    Si f(a)⋅f(c)<0, la raíz está en [a,c], actualiza b=c
    
    Si f(c)⋅f(b)<0f(c)⋅f(b)<0, la raíz está en [c,b], actualiza a=c
    
    
    
    Repite el proceso hasta que el intervalo sea suficientemente pequeño, es decir:
        b -a /2 < tolerancia

"""

########
from sympy import symbols

x = symbols('x')

#intervalo
a = 1
b = 2


#funcion
def f(x):
    return x**3 + 4*x**2 - 10


def biseccionMetodo_N(f, a, b, n):
    if f(a) * f(b) >= 0:
        print("No se cumple Bolzano: no hay raíz en el intervalo")
        return None

    for i in range(n):
        c = (a + b) / 2
        if f(c) == 0:
            print(f"Raíz exacta encontrada en x = {c}")
            return c
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c

    c = (a + b) / 2
    print("Raíz aproximada tras", n, "iteraciones:", c)
    return c


def biseccionMetodo_tolerancia(f, a, b, tolerancia):
    if f(a) * f(b) >= 0:
        print("No se cumple Bolzano: no hay raíz en el intervalo")
        return None

    while (b - a) / 2 > tolerancia:
        c = (a + b) / 2
        if f(c) == 0:
            print(f"Raíz exacta encontrada en x = {c}")
            return c
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c

    c = (a + b) / 2
    print(f"Raíz aproximada con tolerancia {tolerancia}: {c}")
    return c


            
biseccionMetodo_N(f,a,b,3)
biseccionMetodo_tolerancia(f,a,b,0.0001)




###ejerciico 2



"""
    Ejercicio 2 (Ejemplo 2, Tema 4)
    
    Resuelve la ecuación de Kepler
    x − e sen x − b = 0
    para e = 0.5 y b = 0.7 en el intervalo [0, 2], utilizando los métodos de bisección,
    punto fijo, Newton-Raphson y de la secante. Para ello, escribe los procedimientos 
    correspondientes en python de manera análoga al ejercicio anterior.
    Al igual que en el citado ejercicio, representa la función anterior en el intervalo adecuado
    para comprobar si se verifican las hipótesis del teorema de Bolzano.
    Por último, verifica las tablas mostradas en el ejemplo de clase, y compara con el valor 
    exacto que
    proporciona el comando adecuado de python.

"""

# Ejercicio 2 - Resolución de la ecuación de Kepler con métodos numéricos
# x - e*sin(x) - b = 0, con e = 0.5 y b = 0.7 en el intervalo [0, 2]

import math
import matplotlib.pyplot as plt

# Parametros dados
e = 0.5
b = 0.7

# Definimos la funcion f(x) = x - e*sin(x) - b
def f(x):
    return x - e * math.sin(x) - b

# Derivada de f(x), necesaria para Newton-Raphson
def f_prime(x):
    return 1 - e * math.cos(x)

# Funcion para graficar y comprobar Bolzano
def graficar_funcion():
    xs = [i * 0.01 for i in range(0, 201)]  # de 0 a 2
    ys = [f(x) for x in xs]
    plt.plot(xs, ys)
    plt.axhline(0, color='gray', linestyle='--')
    plt.title("Gráfica de f(x) = x - e*sin(x) - b")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid(True)
    plt.show()

# Metodo de Biseccion
def biseccion(f, a, b, tolerancia):
    if f(a) * f(b) >= 0:
        print("No se cumple Bolzano")
        return None

    while (b - a) / 2 > tolerancia:
        c = (a + b) / 2
        if f(c) == 0:
            return c  # Raiz exacta
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2

# Metodo de Punto Fijo: x = g(x)
def punto_fijo(g, x0, tolerancia, max_iter):
    for i in range(max_iter):
        x1 = g(x0)
        if abs(x1 - x0) < tolerancia:
            return x1
        x0 = x1
    return x0

# Metodo de Newton-Raphson
def newton_raphson(f, f_prime, x0, tolerancia, max_iter):
    for i in range(max_iter):
        x1 = x0 - f(x0) / f_prime(x0)
        if abs(x1 - x0) < tolerancia:
            return x1
        x0 = x1
    return x0

# Metodo de la Secante
def secante(f, x0, x1, tolerancia, max_iter):
    for i in range(max_iter):
        if f(x1) - f(x0) == 0:
            return None  # Evitar division por 0
        x2 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0))
        if abs(x2 - x1) < tolerancia:
            return x2
        x0, x1 = x1, x2
    return x1

# Funcion g(x) para punto fijo --> x = e*sin(x) + b
def g(x):
    return e * math.sin(x) + b

# Graficamos para verificar Bolzano
graficar_funcion()

# Aplicamos cada metodo
print("\nMetodo de Biseccion:")
raiz_biseccion = biseccion(f, 0, 2, 0.0001)
print("Raiz aproximada:", raiz_biseccion)

print("\nMetodo de Punto Fijo:")
raiz_pf = punto_fijo(g, 1.0, 0.0001, 100)
print("Raiz aproximada:", raiz_pf)

print("\nMetodo de Newton-Raphson:")
raiz_newton = newton_raphson(f, f_prime, 1.0, 0.0001, 100)
print("Raiz aproximada:", raiz_newton)

print("\nMetodo de la Secante:")
raiz_secante = secante(f, 0.5, 1.5, 0.0001, 100)
print("Raiz aproximada:", raiz_secante)

# Verificamos con solucion exacta usando scipy
from scipy.optimize import fsolve
raiz_exacta = fsolve(lambda x: x - e * math.sin(x) - b, 1)[0]
print("\nValor exacto con fsolve:", raiz_exacta)

#################################################################################
######################sistemas lineales##########################################
"""
Ejercicio (Ejemplo Introducción, Tema 5)

a) programe las subrutinas de factorización, sustitución progresiva y sustitución 
   regresiva que permitan resolver sin usar pivotaje un sistema lineal n × n


b) aplique el apartado anterior para resolver el siguiente sistema lineal



    [2, 3,   4,  5]      [x1]     [ 5]
    [6, 15, 19, 23]  *   [x2]  =  [30]
    [8, 42, 60, 70]     [x3]      [98]
    [12, 60, 1, 17]     [x4]     [144]
    
    
c) tilícese para calcular la matriz inversa de A y su determinante, 
    siendo A la matriz del sistema
    del apartado anterior.
    
"""
# a)
    def factorizacion_LU_sin_pivotaje(A):
    """
    Esta función realiza la factorización LU de una matriz A sin pivotaje.
    Devuelve dos matrices: L (triangular inferior) y U (triangular superior)
    tal que A = L * U.
    """
    n = A.shape[0]  # Obtenemos el número de filas (la matriz debe ser cuadrada)
    L = eye(n)      # L empieza como matriz identidad (1s en la diagonal)
    U = A.copy()    # U empieza como copia de A, y se irá transformando

    for k in range(n):
        pivot = U[k, k]
        if pivot == 0:
            raise ValueError("No se puede continuar sin pivotaje, pivote cero.")
        
        # Recorremos todas las filas por debajo del pivote actual
        for i in range(k + 1, n):
            m = U[i, k] / pivot      # Calculamos el multiplicador para eliminar la columna
            L[i, k] = m              # Lo guardamos en L
            U[i, :] = U[i, :] - m * U[k, :]  # Restamos m * fila_k a la fila_i

    return L, U

#b)
from sympy import Matrix, eye

def sustitucion_progresiva(L, b):
    """
    Resuelve L * y = b cuando L es triangular inferior.
    """
    n = L.shape[0]
    y = [0]*n
    for i in range(n):
        suma = sum(L[i,j]*y[j] for j in range(i))
        y[i] = (b[i] - suma) / L[i,i]
    return Matrix(y)


def sustitucion_regresiva(U, y):
    """
    Resuelve U * x = y cuando U es triangular superior.
    """
    n = U.shape[0]
    x = [0]*n
    for i in reversed(range(n)):
        suma = sum(U[i,j]*x[j] for j in range(i+1, n))
        x[i] = (y[i] - suma) / U[i,i]
    return Matrix(x)




# Matriz del sistema
A = Matrix([
    [2, 3, 4, 5],
    [6, 15, 19, 23],
    [8, 42, 60, 70],
    [12, 60, 1, 17]
])

# Vector del lado derecho
b = Matrix([5, 30, 98, 144])

# Factorización LU sin pivotaje
L, U = factorizacion_LU_sin_pivotaje(A)

# Sustitución progresiva: resolvemos L*y = b
y = sustitucion_progresiva(L, b)

# Sustitución regresiva: resolvemos U*x = y
x = sustitucion_regresiva(U, y)

# Mostramos el resultado
print("Solución del sistema:")
print(x)

#----apartado c; matriz inversa
def determinante_por_LU(U):
    """
    Calcula el determinante como el producto de la diagonal de U.
    (L no cambia el determinante porque su diagonal son todos unos)
    """
    det = 1
    for i in range(U.shape[0]):
        det *= U[i, i]
    return det


def matriz_inversa_por_LU(L, U):
    """
    Calcula la inversa de A resolviendo A * X = I
    usando LU sin pivotaje.
    """
    n = L.shape[0]
    inversa = []
    
    # Creamos matriz identidad
    I = eye(n)
    
    for i in range(n):
        # Tomamos cada columna de la identidad
        e = I[:, i]
        # Resolvemos L * y = e
        y = sustitucion_progresiva(L, e)
        # Luego resolvemos U * x = y
        x = sustitucion_regresiva(U, y)
        # Agregamos columna a la inversa
        inversa.append(x)
        
    # Formamos la matriz inversa como columnas
    return Matrix.hstack(*inversa)


# Usamos las matrices ya obtenidas:
# L, U = factorizacion_LU_sin_pivotaje(A)

# Determinante
det_A = determinante_por_LU(U)
print(f"Determinante de A: {det_A}")

# Inversa
inversa_A = matriz_inversa_por_LU(L, U)
print("Matriz inversa de A:")
print(inversa_A)
