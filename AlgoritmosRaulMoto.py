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
import sympy as sp



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

######################################################################################
###################################RNZ###################################################
######################################################################################
##### EJERCICIO 1

# DECLARAMOS LA FUNCION
def f(x):
    return x**3 + 4*x**2 - 10

## BISECCION CON ITERACION
def biseccion_N(f, a, b, N):
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) y f(b) deben tener signos opuestos.")
        
    print(f"n {' ':<4} an {' ':<15} bn {' ':<15} cn {' ':<15} |an-bn|")
    print("-" * 70)
    
    for n in range(N + 1):
        c = (a + b) / 2
        fc = f(c)
        print(f"{n}{' ':<5} {a:<12.6f} {' ':<5} {b:<12.6f} {' ':<5} {c:<12.6f} {' ':<5} {abs(b-a):<12.6f}")

        if fc == 0:
            return c  # ¡Raíz exacta encontrada!
        elif f(a) * fc < 0:
            b = c
        else:
            a = c

    return c

raiz = biseccion_N(f, 1, 2, 13)
print("Raíz aproximada= " + str(raiz))

def biseccion_TOL(f, a, b, TOL):
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) y f(b) deben tener signos opuestos.")
    i = 0
    continuar = True
    extra = False
    print(f"n {' ':<3} an {' ':<9} bn {' ':<9} cn {' ':<9} f(cn) {' ':<8} |an-bn|")
    print("-" * 69)
    while continuar:
        c = (a + b) / 2
        fc = f(c)
        an_bn = (b - a) / 2
        print(f"{i} {a:12.6f} {' ':<3} {b:<12.6f} {c:<12.6f} {fc:<14.6f} {an_bn:<12.6f}")
        i += 1
        if fc == 0:  # Raíz exacta encontrada
            break
        elif f(a) * fc < 0:
            b = c
        else:
            a = c
        if an_bn <= TOL:
            if not extra:
                extra = True
            else:
                continuar = False
    return c
raiz = biseccion_TOL(f, 1, 2, 0.01)
print("Raíz aproximada= " + str(raiz))

##### EJERCICIO 2
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 01:35:08 2025

@author: renzo
"""

import numpy as np
import matplotlib.pyplot as plt

# Parámetros
e = 0.5
b = 0.7

# Función de Kepler
def f(x):
    return x - e * np.sin(x) - b

# Gráfico en [0, 2]
x_vals = np.linspace(0, 2, 400)
y_vals = f(x_vals)

plt.plot(x_vals, y_vals, label='f(x)')
plt.axhline(0, color='orange', linestyle='-')
plt.title("Ecuación de Kepler: f(x) = x - 0.5·sin(x) - 0.7")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.legend()
plt.show()

## BISECCION CON ITERACION
def biseccion_N(f, a, b, N):
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) y f(b) deben tener signos opuestos.")
        
    print(f"n {' ':<4} an {' ':<15} bn {' ':<15} cn {' ':<15} |an-bn|")
    print("-" * 70)
    
    for n in range(N + 1):
        c = (a + b) / 2
        fc = f(c)
        print(f"{n}{' ':<5} {a:<12.6f} {' ':<5} {b:<12.6f} {' ':<5} {c:<12.6f} {' ':<5} {abs(b-a):<12.6f}")

        if fc == 0:
            return c  # ¡Raíz exacta encontrada!
        elif f(a) * fc < 0:
            b = c
        else:
            a = c

    return c

raiz = biseccion_N(f, 0, 2, 5)
    
print("Raíz aproximada= " + str(raiz))

def biseccion_TOL(f, a, b, TOL):
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) y f(b) deben tener signos opuestos.")
    i = 0
    continuar = True
    extra = False
    print(f"n {' ':<3} an {' ':<9} bn {' ':<9} cn {' ':<9} f(cn) {' ':<8} |an-bn|")
    print("-" * 69)
    while continuar:
        c = (a + b)/ 2
        fc = f(c)
        an_bn = (b - a)/ 2
        print(f"{i} {a:12.6f} {' ':<3} {b:<12.6f} {c:<12.6f} {fc:<14.6f} {an_bn:<12.6f}")
        i+=1
        if fc == 0:
            break
        elif f(a)*fc < 0:
            b=c
        else:
            a=c
        if an_bn <= TOL:
            if not extra:
                extra=True
            else:
                continuar=False
    return c

raiz = biseccion_TOL(f, 0, 2, 0.001)
print("Raíz aproximada= " + str(raiz))

###############################################################################

##PUNTO FIJO
# Función g(x) para punto fijo
def g(x):
    return e * np.sin(x) + b

def punto_fijo(g, x0, N):
    print(f"{'n':<3} {'x_n':<12}")
    print("-"*50)
    guardar_r = 0
    for n in range(0, N+1):
        x1 = g(x0)
        print(f"{n:<3} {x0}")
        guardar_r = x0
        x0 = x1
    return guardar_r

print("---PUNTO FIJO---")
print("Punto fijo raiz= " + str(punto_fijo(g, 0, 5)))

def puntofijo_tol(g, x0, TOL):
    i = 0
    while True:
        x1 = g(x0)
        error = abs(x1 - x0)
        print(f"Iteración {i}: x = {x1}, error = {error}")
        if error < TOL:
            break
        x0 = x1
        i += 1
    return x1
print("\n--- PUNTO FIJO TOL ---")
print("Punto fijo tol raiz= " + str(puntofijo_tol(g, 0, 0.001)))
###############################################################################


def newton_RaphsonN(f, x0, N):
    h = 1e-6  # Paso pequeño para derivada numérica
    print(f"      xn")
    print(f"0      {x0}")
    for i in range(1, N + 1):
        f_x0 = f(x0)
        # Derivada aproximada: f'(x) ≈ (f(x + h) - f(x)) / h
        df_x0 = (f(x0 + h) - f_x0) / h

        if df_x0 == 0:
            raise ZeroDivisionError("La derivada numérica es cero. No se puede continuar.")

        x1 = x0 - f_x0 / df_x0
        print(f"{i}      {x1}")
        x0 = x1

    return x0

print("\n--- Newton-Raphson ---")
print("Newton-Raphson raiz= " + str(newton_RaphsonN(f, 0, 4)))
###############################################################################
##SECANTE
def secanteN(f, x0, x1, N):
    print(f"{' ':<5} xn {' ':<22} |xn-xn-1|")
    print(f"0 {' ':<5} {x0} {' ':<22} 0")
    print(f"1 {' ':<5} {x1} {' ':<22} {abs(x1 - x0)}")

    for i in range(2, N + 1):
        f_x0 = f(x0)
        f_x1 = f(x1)
        if f_x1 - f_x0 == 0:
            raise ZeroDivisionError("División por cero en la fórmula de la secante.")
        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
        error = abs(x2 - x1)
        print(f"{i} {' ':<5} {x2} {' ':<5} {error}")

        # Actualizar valores para la siguiente iteración
        x0, x1 = x1, x2
    return x2
secanteN(f, 0, 2, 4)

def secanteTOL(f, x0, x1, TOL):
    i = 1
    error = abs(x1 - x0)
    print(f"{' ':<7} xn {' ':<22} |xn-xn-1|")
    print(f"0 {' ':<5} {x0} ")
    print(f"1 {' ':<5} {x1} {' ':<22} {error}")

    while error > TOL:
        f_x0 = f(x0)
        f_x1 = f(x1)

        if f_x1 - f_x0 == 0:
            raise ZeroDivisionError("División por cero en la fórmula de la secante.")

        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
        error = abs(x2 - x1)

        i += 1
        print(f"{i} {' ':<5} {x2} {' ':<17} {error}")

        x0, x1 = x1, x2

    return x2
print("\n--- SECANTE TOL---")
secanteTOL(f, 0, 2, 0.001)


######################################################################################
###################################RNZ###################################################
######################################################################################






#--------------------------------------DICCIONARIO---------------------------------------
#--------------------------------------DICCIONARIO---------------------------------------
#--------------------------------------DICCIONARIO---------------------------------------
#--------------------------------------DICCIONARIO---------------------------------------




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

#--------------------------------------DICCIONARIO fin---------------------------------------
#--------------------------------------DICCIONARIO fin---------------------------------------
#--------------------------------------DICCIONARIO fin---------------------------------------
#--------------------------------------DICCIONARIO fin---------------------------------------


#--------------------------------------Tema2 -------------------------------------------
#--------------------------------------Tema2 -------------------------------------------
#--------------------------------------Tema2 -------------------------------------------


# -*- coding: utf-8 -*-
"""

@author: I Gomara based on JI Farran notes
"""

#PRÁCTICA 2: INTERPOLACIÓN POLINÓMICA  

# Vamos a obtener los polinomios de Taylor para interpolar  
# Usamos la librería sympy  

from sympy import *  

#definimos x como incógnita

x = symbols("x") 

# Desarrollo en serie de Taylor hasta orden 3 en x=1  
#función series de numpy (x es incógnita, 
#1 punto donde se construye Taylor, 4 número de términos - orden 3)

s = (x * exp(x)).series(x, 1, 4) 
print(s) 

r= (sin(x)).series(x, 0, 4) 
print(r) 

# Obtenemos una aproximación eliminando O() y sustituyendo
# 0() es el error (resto), hay que quitarlo para evaluar polinomio  

s = s.removeO() 
print(s)

s1=s.subs(x, 3) 
float(s1)

# Evaluación númerica (no simbólica) con evalf  

s.evalf(subs={x: 3}) 

# Se podrían ir calculando los coeficientes de Taylor y construyendo el polinomio  

from sympy import diff, factorial 

f = x * exp(x) 
a0= f.subs(x, 1) 
print(a0)
a1 = diff(f, x).subs(x, 1) 
print(a1)
a2 = diff(f, x, 2).subs(x, 1) / factorial(2) 
print(a2)
a3 = diff(f, x, 3).subs(x, 1) / factorial(3) 
print(a3)
P = a0 + a1*(x-1) + a2*(x-1)**2 + a3*(x-1)**3 
print(P) 
P.subs(x, 3) 
P.evalf(subs={x: 3}) 

P==s

# Interpolación polinómica automática con numpy (lineal) 

import numpy as np 
import sympy as sp
import matplotlib.pyplot as plt

# nodos y valores de la función  

xp = [1, 2, 3] 
fp = [3, 2, 0] 

# Interpolamos en el punto x = 2.5  

np.interp(2.5, xp, fp) 

# Recordar que se pueden dibujar funciones con matplotlib.pyplot  

#esto es para discretizar el eje x (cuanta mayor resolución más suave el ploteado)

t=np.arange(1,3.1,0.1)
t2=np.linspace(1,3,21)

print(t)
print(t2)

#vamos a dibujar funcion seno

ejex=np.linspace(-np.pi, np.pi, 100);

#comando para que funciones de sympy admitan evaluación introduciendo un array 
#en lugar de un dato

fsin = sp.lambdify(x, sp.sin(x), 'numpy')
fcos = sp.lambdify(x, sp.cos(x), 'numpy')

#ahora ya se puede introducir el array en la función, y devuelve un array (las f(y))

ejey1=fsin(ejex)
ejey2=fcos(ejex)

print(ejex)
print(ejey1)

plt.plot(ejex,ejey1,label="sin(x)")
plt.plot(ejex,ejey2,label="cos(x)")
plt.legend()
plt.show()
    
#########################
#ANOTACIONES ADICIONALES SOBRE FUNCIONES - REPASO


#Una función en python es un conjunto de intrucciones escritas en código cuya 
#finalidad es realizar una tarea específica.

#las funciones generalmente poseen uno o varios input y uno o varios output


#EJEMPLO, VAMOS A CONSTRUIR UNA FUNCIÓN QUE REALICE LA SUMA DE 3 NÚMEROS

import numpy as np
import sympy as sp
import math as mt

def suma3(a,b,c):
    """Esta función realiza la suma de los valores a, b y c que indica el usuario"""
    suma=a+b+c
    return [a, b, c, suma]

#Si se ejecuta hasta aquí la función no hace nada. Python lee las instrucciones
#de la misma y la compila. Para ejecutarla, ha de ser llamada, indicando los valores
#de a, b y c (hasta ahora variables genéricas)

suma3(1,2,3)[3]


#si ejecutamos de este modo la función devuelve los output, pero no los guarda en una variable

output=suma3(1,2,3)
output

#ahora ya sí hemos guardado los output.

#SE DESACONSEJA LA COMUNICACIÓN POR CONSOLA CON LA FUNCIÓN UTILIZANDO EL COMANDO INPUT

#supongamos que creamos la misma función, pero indicamos los datos por consola

def suma3_input():
    """Esta función realiza la suma de los valores a, b y c que indica el usuario"""
    a = input("Introduzca el primer número: \n")
    a=float(a)
    b = input("Introduzca el segundo número: \n")
    b=float(b)
    c = input("Introduzca el tercer número: \n")
    c=float(c)
    suma=a+b+c
    return [a, b, c, suma]

suma3_input()

#como vemos, realiza la misma función, pero es mucho menos versátil, ya que se precisa
#comunicación por pantalla.

#supongamos ahora que se pide calcular la suma de 3 números consecutivos, empezando por
#el 1 y terminando en el 100

#con la primera versión es bastante directo

sumas=[]
for i in range(1,101):
    sumas.append(suma3(i, i+1, i+2)[3])
sumas
    
#pero con la segunda sería inviable, habría que meter 300 números a mano

sumas=[]
for i in range(1,101):
    sumas.append(suma3_input()[3])   
    
#POR TANTO, PARA EVITAR PROBLEMAS EN ESTE SENTIDO, SE RECOMIENDA UTILIZAR SIEMPRE
#LA PRIMERA VERSIÓN (A NO SER QUE SE INDIQUE EXPRESAMENTE LO CONTRARIO)

#TRABAJO CON TABLAS

#se puede crear una tabla con 0 o 1s

Y=np.ones((2,2))
Y

T=np.zeros((3,4))

#para meter datos en filas, columnas y posiciones

T

#filas

T[:,0]=[1,2,3]    
T

#filas a partir de una posición hasta el final

T[1,1:]=[2,2,2]

T

#filas desde posición inicial hasta final (recordamos que siempre que queda en 
#la anterior a la final)

T[2,2:4]=[3,3]

T
#--------------------------------------Tema2 -------------------------------------------
#--------------------------------------Tema2 -------------------------------------------
#--------------------------------------Tema2 -------------------------------------------



#--------------------------------------Tema3 -------------------------------------------
#--------------------------------------Tema3 -------------------------------------------
#--------------------------------------Tema3 -------------------------------------------

# -*- coding: utf-8 -*-
"""

@author: i_gom based on JI Farrán notes
"""

#PRÁCTICA 3: DERIVACIÓN E INTEGRACIÓN NUMÉRICAS  

# algo de Álgebra Lineal   

# usaremos la librería numpy  

import numpy as np  

# para trabajar con matrices en lugar de arrays multidimensionales

A = np.zeros((4, 4)) 
A


A1=np.zeros((4,4))
A1

# otras opciones; ones, eye, rand ...  

A2=np.eye(4,4)
A2

A[0, 0] = 1 
print(A) 

# ojo además de copiar la matriz se copian las operaciones que se hagan en ella  

B = A 

# si solo queremos copiar y que se puedan hacer operaciones independientemente  

C = np.copy(A) 

A[0, 0] = 3

A
B
C

# cálculo matricial  

#crear array bidimensional
X = np.array([[1, 2], [3, 4]]) 
X

#convertir a matriz (algunas diferencias en funcionalidades)
M = np.asmatrix(X) 
M

#crear matriz
N = np.matrix([[1, 1], [0, 1]]) 
N

#dimensiones de M (filas, columnas)
M.shape 

#productos escalares/matriciales
np.dot(2, 3) 

np.dot(M, 2) 

np.dot(2, M) 

np.dot(M,M)
np.dot(X, X) 

#otra opción
np.matmul(M, M) 

#ojo que * es producto matricial si M es matriz
M * M 

M + M 

M - M 

3 * M 

M * 3 

#traspuesta
np.transpose(M) 

M
M.T 

#traza

A
np.trace(A) 

#diagonal principal

np.diag(A) 

# producto escalar  

a = np.array([1,2,3]) 
b = np.array([0,1,0]) 
np.inner(a, b) 
np.dot(a,b) #lo mismo si tenemos dos vectores

# concatenar por filas (mismo número de columnas)  

X
np.concatenate((X, X)) 
np.concatenate((X, X),axis=0) 

# concatenar por columnas (mismo número de filas)  
#recordemos que axis 0 es filas y 1 columnas

np.concatenate((X, X), axis=1) 

#SISTEMAS LINEALES

from numpy import linalg as LA 

print(M) 

#calcula determinante
LA.det(M) 

#calcula inversa
LA.inv(M) 

#calcula rango
LA.matrix_rank(M) 

# resolución de sistemas lineales  

a = np.array([[3,1], [1,2]]) 
b = np.array([9,8]) 

x = LA.solve(a, b) 
print(x) 

# ejemplo: generar una matriz de tipo Vandermonde  

def vandermonde(nodos): 
    """requiere haber importado numpy as np 
    #nodos es un array con los nodos"""
    n = len(nodos) 
    X = np.zeros((n, n)) 
    for j in range(n): #iteración columnas
        for i in range(n): #iteración filas
            X[i, j] = nodos[j] ** i #coge el punto xj y lo eleva al indice fila
    M = np.asmatrix(X) #convierte a matrix
    return(M) #devuelve matrix

vandermonde([2,3,4])

# la función anterior se usará para fórmulas de cuadratura y derivación numéricas  
# mediante el método de coeficientes indeterminados  

# ejemplo: aproxima la derivada del sen(x) en x0=0  

#con una fórmula central de 5 puntos para h=0.1

h = 0.1 #distancia
x0 = 0. #la coma al final para que sea float
P = [x0 - 2*h, x0 - h, x0, x0 + h, x0 + 2*h]  #nodos
A = vandermonde(P) #calculamos vadermonde para esos nodos (recordar diapositivas clase)
print(A)
b = np.array([0, 1, 2*x0, 3*x0**2, 4*x0**3]) #vector al otro lado de la igualdad (clase)

alphas = LA.solve(A, b) #resolvemos el sistema para obtener los alfa-sub-i
print(alphas) #valores de los coeficientes al resolver el sistema

#ahora hay que multiplicar los coeficientes (alphas) por el valor de los f(xi)
#recordemos D5(x)=alphas(i)*f(xi)

import sympy as sp
x = sp.symbols("x") 

F = 0. 
for i in range(5): 
    F = F + alphas[i] * sp.sin(P[i]) 
print(F) 

#otra opción para F, usar el producto vectorial

fsin = sp.lambdify(x, sp.sin(x), 'numpy')
fxi=fsin(P) #P son los 5 nodos, calculamos los f(xi)

F1=np.dot(alphas,fxi)
print(F1)

F==F1

# cuadratura numérica automática con python  

from scipy import integrate 
import sympy as sp
x=sp.symbols("x")

#Fundamental algorithms for scientific computing in Python https://scipy.org/

x2 = lambda x: x**2 

#recordemos el lambdify
#his module provides convenient functions to transform SymPy expressions to lambda functions 
#which can be used to calculate numerical values very fast

x2p = sp.lambdify(x, x**2, 'numpy')

# devuelve una aproximación numérica y el orden del error  

integrate.quad(x2, 0, 4) 
integrate.quad(x2p, 0, 4) 

# comparamos con el valor exacto (integral de x2 -> x3/3)

print(4**3 / 3.) 

# valor exacto de una integral con sympy (x2)  

import sympy as sp
F = sp.integrate(x**2, x)  # paciencia, lleva unos segundos 
F

#probamos con cualquier otra función

F2 = sp.integrate(sp.sin(x), x) 
F2

F3 = sp.integrate(x**5, x) 
F3

#sustitución en un intervalo-regla de Barrow 

F3.subs(x, 2) - F3.subs(x, 0)  # valor exacto

#también se puede con integrate especificando intervalo (más rápido)

F3 = sp.integrate(x**5, (x, 0, 2))
F3
float(F3)  

F.evalf(subs={x: 2.}) - F.evalf(subs={x: 0.})  # aproximación numérica  

################################NOTACIONES#################################

import sympy as sp

def integral_definida(expr, variable, a, b):
    """
    Calcula la integral definida de una función entre dos límites.

    Parámetros:
    - expr: la expresión simbólica (por ejemplo, x**2)
    - variable: la variable de integración (por ejemplo, x)
    - a: límite inferior (por ejemplo, 0)
    - b: límite superior (por ejemplo, 2)

    Ejemplo de uso:
    x = sp.symbols('x')
    integral_definida(x**2, x, 0, 2)
    """
    # Calculamos la integral definida desde a hasta b
    resultado = sp.integrate(expr, (variable, a, b))

    # Mostramos paso a paso
    print(f"Integral de {expr} desde {a} hasta {b} respecto a {variable}: {resultado.evalf()}")

    # Retornamos el resultado numérico
    return resultado.evalf()
#################################################################

def integral_indefinida_en_punto(expr, variable, punto):
    """
    Calcula la integral indefinida de una función y evalúa su valor en un punto dado.

    Parámetros:
    - expr: la expresión simbólica (por ejemplo, sin(x))
    - variable: la variable de integración (por ejemplo, x)
    - punto: el valor numérico donde se quiere evaluar (por ejemplo, pi)

    Ejemplo de uso:
    x = sp.symbols('x')
    integral_indefinida_en_punto(sp.sin(x), x, sp.pi)
    """
    # Calculamos la integral indefinida (sin límites)
    integral = sp.integrate(expr, variable)

    # Evaluamos la integral en el punto dado
    valor_en_punto = integral.subs(variable, punto)

    # Mostramos los resultados paso a paso
    print(f"Integral indefinida de {expr} respecto a {variable}: {integral}")
    print(f"Evaluada en el punto {punto}: {valor_en_punto.evalf()}")

    # Retornamos el valor evaluado numérico
    return valor_en_punto.evalf()


# Definimos la variable simbólica
x = sp.symbols('x')

# 1. Integral definida de x^2 de 0 a 2
resultado1 = integral_definida(x**2, x, 0, 2)

# 2. Integral indefinida de sin(x), evaluada en pi
resultado2 = integral_indefinida_en_punto(sp.sin(x), x, sp.pi)




import sympy as sp

def taylor(f, x, a, n):
    """
    Calcula el polinomio de Taylor de una función f(x), centrado en el punto a, con n términos.

    Parámetros:
    - f: función simbólica. Ejemplo: sp.exp(x), sp.sin(x), x**2 + x
    - x: variable simbólica. Ejemplo: x = sp.symbols('x')
    - a: punto alrededor del cual se construye la serie de Taylor. Ej: a = 0
    - n: número de términos del polinomio. Ej: 4 → calcular hasta la 3ª derivada

    Ejemplo de uso:
    x = sp.symbols('x')
    taylor(sp.exp(x), x, 0, 4)
    """

    # Inicializamos el polinomio en 0
    taylor_pol = 0

    print(f"\nPolinomio de Taylor de f(x) = {f} centrado en x = {a} hasta {n} términos:\n")

    for i in range(n):
        # Calculamos la i-ésima derivada de f
        derivada_i = sp.diff(f, x, i)

        # Evaluamos la derivada en el punto a
        derivada_evaluada = derivada_i.subs(x, a)

        # Calculamos el término de la serie de Taylor: f^(i)(a)/i! * (x - a)^i
        termino = (derivada_evaluada / sp.factorial(i)) * (x - a)**i

        # Lo sumamos al polinomio
        taylor_pol += termino

        # Mostramos el paso
        print(f"Término {i}: ({derivada_evaluada})/{i}! * (x - {a})^{i} = {termino}")

    # Mostramos el resultado final
    print(f"\nPolinomio de Taylor completo: {sp.simplify(taylor_pol)}")

    return taylor_pol



# Definimos la variable simbólica
x = sp.symbols('x')

# Calculamos el polinomio de Taylor de e^x centrado en 0 con 4 términos
taylor(sp.exp(x), x, 0, 4)

####################################################################

################################################################
#-----------------------lubrerias------------------
from sympy import symbols, sympify, expand, simplify, factor

x = symbols('x')  # Crea una variable simbólica x
expr = sympify("x**2 + 3*x + 2")  # Convierte un string a una expresión simbólica

expand((x + 1)**2)       # x**2 + 2*x + 1
simplify(x**2 + 2*x + 1) # x**2 + 2*x + 1 -> intenta reducirla
factor(x**2 + 2*x + 1)   # (x + 1)**2


from sympy import diff, integrate

diff(expr, x)                 # Derivada de expr respecto a x
integrate(expr, x)           # Integral indefinida respecto a x
integrate(expr, (x, 0, 2))   # Integral definida de x=0 a x=2


from sympy import sin, cos, tan, exp, log, sqrt, factorial

sin(x), cos(x), tan(x)  # Trigonométricas
exp(x)                  # e^x
log(x)                  # logaritmo natural (base e)
log(x, 10)              # logaritmo en base 10
sqrt(x)                 # raíz cuadrada
factorial(5)            # 5! = 120

#------------------
from sympy import series, summation, Symbol

series(exp(x), x, 0, 5)  # Taylor de e^x en x=0 hasta término x^4
i = Symbol('i')
summation(i**2, (i, 1, 3))  # 1^2 + 2^2 + 3^2 = 14


expr.subs(x, 2)    # Sustituye x por 2
expr.evalf()       # Evalúa a número decimal
from sympy import N
N(expr)            # Igual que evalf()



from sympy import Matrix

M = Matrix([[1, 2], [3, 4]])
M.det()        # Determinante
M.inv()        # Inversa
M.rref()       # Forma reducida por filas



# Derivar
diff(f, x)

# Integrar
integrate(f, x)  # Indefinida
integrate(f, (x, a, b))  # Definida

# Sustituir valores
f.subs(x, valor)

# Evaluar a decimal
f.evalf()


# Imprime del 0 al 4
for i in range(5):
    print(i)


# Imprime del 2 al 6
for i in range(2, 7):
    print(i)


# Imprime del 0 al 10 de 2 en 2
for i in range(0, 11, 2):
    print(i)


# Imprime del 5 al 1
for i in range(5, 0, -1):
    print(i)


frutas = ["manzana", "pera", "uva"]
for i, fruta in enumerate(frutas):
    print(f"{i}: {fruta}")

#enumerate te da el índice y el valor a la vez.

#---------------------------
#-----------bolzano ittt-----
import sympy as sp

def bolzano_biseccion(f, a, b, tol):
    """
    Método de Bisección usando el Teorema de Bolzano
    f: expresión simbólica (ejemplo: x**2 - 4)
    a, b: intervalo inicial donde f(a) * f(b) < 0
    tol: tolerancia (cuándo parar)
    """
    x = sp.symbols('x')
    
    # Evaluamos f(a) y f(b)
    fa = f.subs(x, a)
    fb = f.subs(x, b)

    # Comprobamos que hay cambio de signo
    if fa * fb > 0:
        print("No se puede aplicar Bolzano: f(a) y f(b) tienen el mismo signo.")
        return None

    print("Iteraciones del método de Bolzano (bisección):")
    
    while abs(b - a) > tol:
        # Punto medio
        c = (a + b) / 2
        fc = f.subs(x, c)

        print(f"a = {a}, b = {b}, c = {c}, f(c) = {fc}")

        # Si la raíz está en [a, c]
        if fa * fc < 0:
            b = c
            fb = fc
        # Si la raíz está en [c, b]
        else:
            a = c
            fa = fc

    # Valor final aproximado
    raiz = (a + b) / 2
    print(f"\nAproximación final de la raíz: {raiz}")
    return raiz

x = sp.symbols('x')
f = x**2 - 4  # Esta función tiene una raíz en x=2
bolzano_biseccion(f, 0, 3, 0.001)
