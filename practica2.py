# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""
    Ejercicio 1.
    
    Mediante un bucle, haz una lista con todos los números enteros entre 1 y 100 que sean
    múltiplos de 5 y cuyo cuadrado esté comprendido entre 200 y 5000. Dentro de este bucle escribe una
    variable que cuente los elementos de dicha lista, y comprueba que este contador coincide con la longitud
    de la lista generada

"""
lista = []
count = 0
i = 0
def resolver():
    for i in range(100):
        if i%5 == 0 and i**2 in range(200,5000):
            lista.append(i)
            count = len(lista)
    return count


"""
    Calculense los polinomios de Taylor para las siguientes funciones
    (a) sen x, de grado 2n + 1 en x0 = 0

"""

from sympy import * 
x = symbols("x")
##lista de funciones

## serie: el polinomio  o serie de Tylor
def serieDeTylor(f, x0, n):
    s = f.series(x, x0, n+1).removeO() 
    return s

serieDeTylor(sin(x),0,5)


"""
Ejercicio 1.

    Usando la librería sympy de python, calcula los tres primeros polinomios de Taylor (hasta
    grado 5) para la función f (x) = sen(x) en x = 0. Haz una gráfica conjunta de los tres polinomios y la
    función original para comprobar visualmente el grado de aproximación de cada uno d

"""
import matplotlib.pyplot as plt

def serieDeTylor(f, x0, n):
    s = f.series(x, x0, n+1).removeO() 
    return s

serieDeTylor(sin(x),0,3)#fun,punto,rango

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

"""

Ejercicio 2 (Ejemplo 6, tema 2) Escribe un procedimiento en python para calcular el polinomio inter-
polador de orden n en la forma de Lagrange para una función f dada. Aplíquese a la función f (x) = 1/x
en los nodos.

x0 = 2, x1 = 5/2, x2 = 4

y aproxime el valor de f(3)
"""
import numpy as np
x = symbols("x")

def polinomioInterpolador(funcion):
    xp = [2, 5/2, 4]
    np.interp(2.5, xp, funcion) 
    print(np.interp(2.5, xp, funcion) )
    
funcion = 1/x
polinomioInterpolador(funcion)
