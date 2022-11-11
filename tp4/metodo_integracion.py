from array import array
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np            
import sympy as sym
from sympy.ntheory.continued_fraction import continued_fraction
from sympy.utilities.lambdify import lambdify

#m*h
#m: multiplicador de h
def generarCoeficientesExpansion(cant, m):
    expansion = []
    for i in range(0,cant,1):
        denominador = sym.factorial(i)
        numerador = m**i
        expansion.append(numerador/denominador)
    return np.array(expansion)      

def getDiferenciaExpansionExacta(cant, m_limit_inf, m_limit_sup):

    expansion_exacta_lim_sup = generarCoeficientesExpansion(cant, m_limit_sup) 
    expansion_exacta_lim_inf = generarCoeficientesExpansion(cant, m_limit_inf)

    expansion_exacta_diferencia = expansion_exacta_lim_sup - expansion_exacta_lim_inf
    expansion_exacta_diferencia = np.delete(expansion_exacta_diferencia,0)
    return expansion_exacta_diferencia

def getFormulaAproximacion(m_limit_inf, m_limit_sup):
    exp_h = generarCoeficientesExpansion(10, 1)
    exp_2h = generarCoeficientesExpansion(10, 2)
    exp_exacta = getDiferenciaExpansionExacta(11,m_limit_inf,m_limit_sup)
    
    '''
    v1 = [i1 i2 ...]
    v2 = [e1 e2 ...]
    diff = [f1 f2 ...]

    (i1+a)*b + e1*d = f1
    i2*b + e2*d = f2
    i3*b + e3*d = f3
    '''

    A =np.array([ [float(exp_h[1]), float(exp_2h[1])],
                  [float(exp_h[2]), float(exp_2h[2])] ])
    b = np.array([ float(exp_exacta[1]), float(exp_exacta[2]) ])
    x = np.linalg.solve(A,b)

    a = ((exp_exacta[0] - x[1]*exp_2h[0]) / x[0]) - exp_h[0]

    exp_h[0]+=a
    exp_h *= x[0]
    exp_2h *=x[1]
    print("Coeficientes Expansion exacta:",exp_exacta)
    print("Coeficientes Expansion aproximada:",exp_h+exp_2h)

    if a != 0:
        print("h * ( "+str(x[0])+"*(f(x0+h) + "+str(a)+"*f(x0)) + "+str(x[1])+"*f(x0+2h) )\n")
    else:
        print("h * ( "+str(x[0])+"*f(x0+h) + " +str(x[1])+"*f(x0+2h) )\n")

    


getFormulaAproximacion(1,2.5)
getFormulaAproximacion(0.5,1.5)

