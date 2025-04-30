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
import sympy
from sympy import symbol as sb

x = sb('x');

#intervalo
a = 1
b = 2


#funcion
def f(x):
    res = (x**3 + 4*x**2 - 10)
    return res;

# esta solo iterra N veces
def biseccionMetodo_N(f,a,b,n):
    
    #punto medio
    c = a + b/ 2
    
    #verificamos que se cumple blzano
    if f(a) * f(b) < 0:
        print("Se cumpleBolzano");
    
    for i in range(n):
        if(f(a) * f(c) < 0):
            #la raiz estáentre [a,c]
            b = c
        elif f(c)* f(b) < 0:
            #la raiz estáentre [c,b]
            a = c
    print("raiz aproximada: ",c)
        

# aplicando con nivel de tolerancia
def biseccionMetodo_tolerancia(f,a,b,tolerancia):
    
    #punto medio
    c = a + b/ 2
    
    #verificamos que se cumple blzano
    if f(a) * f(b) < 0:
        print("Se cumpleBolzano");
    
    #biseccion::::
        
    #evaluamos la funcion en c
    while b - a / 2  < tolerancia:
        if(f(a) * f(c) < 0):
            #la raiz estáentre [a,c]
            b = c
        elif f(c)* f(b) < 0:
            #la raiz estáentre [c,b]
            a = c
    print("raiz aproximada: ",c)


            
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
