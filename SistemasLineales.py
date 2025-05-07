# -*- coding: utf-8 -*-
"""
Created on Wed May  7 13:41:24 2025

@author: uva
"""

"""
Ejercicio (Ejemplo Introducción, Tema 5)

a) programe las subrutinas de factorización, sustitución progresiva y sustitución 
   regresiva que permitan resolver sin usar pivotaje un sistema lineal n × n


b) plique el apartado anterior para resolver el siguiente sistema lineal



    [2, 3,   4,  5]      [x1]     [ 5]
    [6, 15, 19, 23]  *   [x2]  =  [30]
    [8, 42, 60, 70]     [x3]      [98]
    [12, 60, 1, 17]     [x4]     [144]
    
    
c) tilícese para calcular la matriz inversa de A y su determinante, 
    siendo A la matriz del sistema
    del apartado anterior.
    
"""

from sympy.matrices import *
#La sustitución progresiva es un método para resolver un sistema de ecuaciones 
#lineales cuando la matriz es triangular inferior


    

def sistemaLineal(M):
    
    #obtenemos priemro el orden de la matriz
    orden = M.shape
    
    #comprobamos si es cuadrada o no
    if orden[0]!=orden[1]:
        print(f"debe ser cuadrada")
        return()
    
    #definimos matriz L, U
    L = eye(orden)
    U = copy(M)
    
    #iterramos sobre pivotes
    for i in range(orden):
        if U[k,k] == 0:
            print(f" no se puede continuar")
            return()
        else:
            for i in range(k +1, orden):
                L[i,k]= U[i,k]/U[k,k]
                U[i,k]= 0
    
    print(f"{U}")

sistemaLineal([[2,2],[2,2]])
                
            
                
        
        
    
    
