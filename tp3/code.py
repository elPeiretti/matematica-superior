from array import array
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np            
import sympy as sym
from sympy.ntheory.continued_fraction import continued_fraction
from sympy.utilities.lambdify import lambdify

data = pd.read_csv(os.getcwd()+"\\cardio_03.csv",",",header=None)

data[0].pop(0)
data[1].pop(0)

ecd_x0 = data[0].tolist()
ecd_y0 = data[1].tolist()

ecd_x = [float(x.replace(',', '.')) for x in ecd_x0]
ecd_y = [float(x.replace(',', '.')) for x in ecd_y0]
#plt.suptitle('SeÃ±al obtenida del ECG')
#plt.grid(True)
#plt.plot(ecd_x,ecd_y)
#plt.show()

'''
f1 = a*x**3 + b*x**2 + c*x + d / 0 - 12.5 - Cuadratica
f2 = e*x**3 + f*x**2 + g*x + h / 12.5 - 30 - Cubica

Condiciones a plantear: Nos quitan grados de libertad
f1(0) = 0 -> d = 0
f1(12.5) - f2(12.5) = 0
f2(30) = 0 
f1'(12.5) - f2'(12.5) = 0

Una vez hecho eso, aplicamos minimos cuadrados para despejar las ecuaciones
(Sumatoria(f1-yreal)^2)' =  0  

Tratamos de plantearlo como una matriz:
[[0, 0, 0, 1, 0, 0, 0, 0],
[12.5**3, 12.5**2, 12.5, 1, -12.5*3, -(12.5)**2, -12.5, -1],
[0, 0, 0, 0, 30**3, 30**2, 30, 1],
[3*(12.5**2), 2*12.5, 1, 0, -3*(12.5**2), -2*12.5, -1, 0 ]]*[a, b, c, d, e, f, g, h]=[0, 0, 0, 0]
'''

# Seria correcto aplicar minimos cuadrados por tramos condicinando los puntos de conexion de las aproximaciones.
# Debemos realizar estas limitaciones a mano y luego aplicar minimos cuadrados con un algoritmo?
# Existe alguna manera de programar estas limitaciones?

#30 - 39 - Cuadratica


def metodo_newton(x, y):
    
    i = 0
    j = 1
    ans = []
    ans.append(y)
    
    for v in range(0,len(x)-1):
        ar = ans[v]
        aux = []
        i = 0
        j = v+1
        while(j < len(x)):
            aux.append((ar[i+1]-ar[i])/(x[j]-x[i]))
            j+=1
            i+=1
        ans.append(aux)
    print(ans)

def metodo_newton2(x, y, dx, dy):
    
    i = 0
    j = 1
    ans = []
    ans.append(y)
    
    aux = []
    i = 0
    j = 1
    while(j < len(x)):
        if(x[j]!=x[i]): 
            aux.append((y[j]-y[i])/(x[j]-x[i]))
        else: 
            aux.append(dy[dx.index(x[i])])
        j+=1
        i+=1
    ans.append(aux)

    for v in range(1,len(x)-1):
        ar = ans[v]
        aux = []
        i = 0
        j = v+1
        while(j < len(x)):
            aux.append((ar[i+1]-ar[i])/(x[j]-x[i]))
            j+=1
            i+=1
        ans.append(aux)
    
    expresion = ''
    contador = 0
    
    for v in ans:
        coef = v[0]
        expresion += str(coef)
        i = 0
        while(i < contador):
            expresion+='(x-' + str(x[i]) + ')'
            i+=1
        expresion+='+'
        contador += 1

    print(expresion)

def metodo_newton3(x, y, dx, dy, d2x, d2y):
    
    i = 0
    j = 1
    ans = []
    ans.append(y)
    
    aux = []
    i = 0
    j = 1
    while(j < len(x)):
        if(x[j]!=x[i]): 
            aux.append((y[j]-y[i])/(x[j]-x[i]))
        else: 
            aux.append(dy[dx.index(x[i])])
        j+=1
        i+=1
    ans.append(aux)

    i=0
    j=2
    aux=[]
    while(j < len(x)):
        if(x[j]!=x[i]): 
            aux.append((ans[1][i+1]-ans[1][i])/(x[j]-x[i]))
        else: 
            aux.append(d2y[d2x.index(x[i])]/2)
        j+=1
        i+=1
    ans.append(aux)

    for v in range(2,len(x)-1):
        ar = ans[v]
        aux = []
        i = 0
        j = v+1
        while(j < len(x)):
            aux.append((ar[i+1]-ar[i])/(x[j]-x[i]))
            j+=1
            i+=1
        ans.append(aux)
    
    expresion = ''
    contador = 0
    
    for v in ans:
        coef = v[0]
        expresion += str(coef)
        i = 0
        while(i < contador):
            expresion+='(x-' + str(x[i]) + ')'
            i+=1
        expresion+='+'
        contador += 1

    #print(expresion)
    z = []
    for i in ans:
       z.append(i[0])
    print(z)

#metodo_newton3([1,2,3,4,4,4,5],[2,3,5,7,7,7,9],[4],[9],[4],[-1])

def todoAFloat(x,y,d1,d2,d3):
    for a in [x,y]:
        for i in range(0,len(a)):
            a[i] = float(a[i])

    for a in [d1,d2,d3]:
        for i in range(0,len(a[0])):
            a[0][i] = float(a[0][i])
        for i in range(0,len(a[1])):
            a[1][i] = float(a[1][i])

    return x,y,d1,d2,d3

def imprimirPolinomio(ans,x):
    expresion = ''
    contador = 0
    
    for v in ans:
        coef = v[0]
        expresion += str(coef)
        i = 0
        while(i < contador):
            expresion+='(x-' + str(x[i]) + ')'
            i+=1
        expresion+='+'
        contador += 1

    print(expresion[:len(expresion)-1]) #le quita el + extra que se agrega.

def metodo_newton_final(x,y,deriv1,deriv2,deriv3):

    x,y,deriv1,deriv2,deriv3 = todoAFloat(x,y,deriv1,deriv2,deriv3)
    
    for d in [deriv1, deriv2, deriv3]:
        for i in range(0,len(d[0])):
            pos = x.index(d[0][i])
            x.insert(pos,d[0][i])
            y.insert(pos,y[pos])

    ans = [y]
    for v in range(0,len(x)-1):
        ar = ans[v]
        aux = []
        i = 0
        j = v+1
        while(j < len(x)):
            if x[j]!=x[i]:
                aux.append((ar[i+1]-ar[i])/(x[j]-x[i]))
            else:                             #hay que setear los valores dados para las derivadas
                if v==0: 
                    aux.append(deriv1[1][deriv1[0].index(x[i])])
                if v==1:
                    aux.append(deriv2[1][deriv2[0].index(x[i])]/2)
                if v==2:
                    aux.append(deriv3[1][deriv3[0].index(x[i])]/6)
            j+=1
            i+=1
        ans.append(aux)

    imprimirPolinomio(ans,x)

#metodo_newton_final([1,2,3,4,5],[2,4,8,16,32], [[1,3,5],[3,5,-20]] ,[[1,3],[-1,70]],[[1,3],[50,2]])
''' los argumentos son:
-el vector de coordenadas x de los puntos por donde tiene que pasar el polinomio [x0,x1,x2,..,xn]
-el vector de coordenadas y de los puntos por donde tiene que pasar el polinomio [y0,y1,y2,...yn]
-un vector con dos vectores, el primero contiene las coordenadas x de la primera derivada y el segundo los valores de las derivadas en las x dadas 
[[x1,x3,x8],[13,4,1]] significa que f'(x1)=13, f'(x3)=4 y f'(x8)=1
-otro vector con dos vectores, el primero contiene las coordenadas x de la segunda derivada y el segundo los valores de las derivadas en las x dadas
- otro un vector con dos vectores, el primero contiene las coordenadas x de la tercera derivada y el segundo los valores de las derivadas en las x dadas

si no se quiere condicionar alguna derivada, directamente hay que pasarle los vectores vacios [[],[]]
'''



'''
f1 = a*x**3 + b*x**2 + c*x + d / 0 - 12.5 - Cuadratica
f2 = e*x**3 + f*x**2 + g*x + h / 12.5 - 30 - Cubica

Condiciones a plantear: Nos quitan grados de libertad
f1(0) = 0 -> d = 0
f1(12.5) - f2(12.5) = 0
f2(30) = 0 
f1'(12.5) - f2'(12.5) = 0

Una vez hecho eso, aplicamos minimos cuadrados para despejar las ecuaciones
(Sumatoria(f1-yreal)^2)' =  0  

Tratamos de plantearlo como una matriz:
[[0, 0, 0, 1, 0, 0, 0, 0],
[12.5**3, 12.5**2, 12.5, 1, -12.5*3, -(12.5)**2, -12.5, -1],
[0, 0, 0, 0, 30**3, 30**2, 30, 1],
[3*(12.5**2), 2*12.5, 1, 0, -3*(12.5**2), -2*12.5, -1, 0 ]]*[a, b, c, d, e, f, g, h]=[0, 0, 0, 0]

[0 - 12.5], recta f1=0
[12.5 - 24.5], parabola f2[a2,b2,c2]
[24.5 - 30], cubica f3[a3,b3,c3,d3]
[30 - 39], parabola f4[a4,b4,c4]
[39 - 57.75], cubica f5[a5,b5,c5,d5]
[57.75 - 66.5], parabola f6[a6,b6,c6]
[66.5 - 90], cubica f7[a7,b7,c7,d7]
[90 - fin], f8=0
'''
def get_condiciones_iniciales(f1,f2,p):
    x = sym.symbols('x')

    if f1 != 0 and f2!=0:
        c1 = f1.subs(x,p) -  f2.subs(x,p)
        c2 = sym.diff(f1,x).subs(x,p) -  sym.diff(f2,x).subs(x,p)

    if f1 == 0:
        c1 = f2.subs(x,p)
        c2 = sym.diff(f2,x).subs(x,p)

    if f2 == 0:
        c1 = f1.subs(x,p)
        c2 = sym.diff(f1,x).subs(x,p)

    return [c1, c2]

def dividirEnIntervalos(ecd_x, ecd_y, puntos):
    intervalos_x = []
    intervalos_y = []
    i=0
    for p in puntos:
        aux_x = []
        aux_y = []
        while p>ecd_x[i]:
            aux_x.append(ecd_x[i])
            aux_y.append(ecd_y[i])
            i+=1
        intervalos_x.append(aux_x)
        intervalos_y.append(aux_y)
    aux_x = []
    aux_y = []
    while(i<len(ecd_x)):
        aux_x.append(ecd_x[i])
        aux_y.append(ecd_y[i])
        i+=1
    intervalos_x.append(aux_x)
    intervalos_y.append(aux_y)

    return intervalos_x,intervalos_y

def calcularSistemaMinimosCuadrados1(coef,phi,ecd_x,ecd_y,puntos):
    x = sym.symbols('x')
    A=[]
    b=[]

    for c in coef:
        fila = []
        for c2 in coef:
            sum=0.0
            for xi in ecd_x:
                try:
                    if xi>=puntos[0] and xi<puntos[1]:
                        sum+=phi[0][c].subs(x,xi)*phi[0][c2].subs(x,xi)
                    if xi>=puntos[1] and xi<puntos[2]:
                        sum+=phi[1][c].subs(x,xi)*phi[1][c2].subs(x,xi)
                    if xi>=puntos[2] and xi<puntos[3]:
                        sum+=phi[2][c].subs(x,xi)*phi[2][c2].subs(x,xi)
                    if xi>=puntos[3] and xi<puntos[4]:
                        sum+=phi[3][c].subs(x,xi)*phi[3][c2].subs(x,xi)
                    if xi>=puntos[4]and xi<puntos[5]:
                        sum+=phi[4][c].subs(x,xi)*phi[4][c2].subs(x,xi)
                    if xi>=puntos[5]and xi<puntos[6]:
                        sum+=phi[5][c].subs(x,xi)*phi[5][c2].subs(x,xi)
                except KeyError:
                    pass

            fila.append(float(sum))
        sum=0.0
        for i in range(len(ecd_x)):
            try:
                if ecd_x[i]>=puntos[0] and ecd_x[i]<puntos[1]:
                    sum+=phi[0][c].subs(x,ecd_x[i])*ecd_y[i]
                if ecd_x[i]>=puntos[1] and ecd_x[i]<puntos[2]:
                    sum+=phi[1][c].subs(x,ecd_x[i])*ecd_y[i]
                if ecd_x[i]>=puntos[2] and ecd_x[i]<puntos[3]:
                    sum+=phi[2][c].subs(x,ecd_x[i])*ecd_y[i]
                if ecd_x[i]>=puntos[3] and ecd_x[i]<puntos[4]:
                    sum+=phi[3][c].subs(x,ecd_x[i])*ecd_y[i]
                if ecd_x[i]>=puntos[4]and ecd_x[i]<puntos[5]:
                    sum+=phi[4][c].subs(x,ecd_x[i])*ecd_y[i]
                if ecd_x[i]>=puntos[5]and ecd_x[i]<puntos[6]:
                    sum+=phi[5][c].subs(x,ecd_x[i])*ecd_y[i]
            except KeyError:
                pass
        b.append(float(sum))
        A.append(fila)
    
    return A,b

def calcularSistemaMinimosCuadrados2(coef,phi,intervalos_x,intervalos_y):
    x = sym.symbols('x')
    A=[]
    b=[]

    for c in coef:
        fila = []
        for c2 in coef:
            sum=0.0
            for i in range(len(intervalos_x)):
                for j in range(len(intervalos_x[i])):
                    try:
                        sum+=phi[i][c].subs(x,intervalos_x[i][j])*phi[i][c2].subs(x,intervalos_x[i][j])
                    except KeyError:
                        pass
            fila.append(float(sum))
            
        sum=0.0
        for i in range(len(intervalos_x)):
                for j in range(len(intervalos_x[i])):
                    try:
                        sum+=phi[i][c].subs(x,intervalos_x[i][j])*intervalos_y[i][j]
                    except KeyError:
                        pass

        b.append(float(sum))
        A.append(fila)
    
    return A,b

def graficarFuncion(tramos,puntos, intervalos_x):
    x = sym.symbols('x')
    intervalos_x[0].append(puntos[0])
    for i in range(1,len(intervalos_x)-1):
        intervalos_x[i].append(puntos[i])

    for i in range(len(tramos)):
        lam_x = sym.lambdify(x,tramos[i],modules=['numpy'])
        x_vals = np.array(intervalos_x[i])
        y_vals = lam_x(x_vals)

        if i==0 or i==len(tramos)-1:
            y_vals = [0 for j in range(len(x_vals))]
        plt.plot(x_vals,y_vals,label='tramo '+str(i),linewidth=2)

    plt.plot(ecd_x,ecd_y,label='muestreada',linewidth=0.5)
    plt.legend(loc='best')
    plt.show()

def minimos_cuadrados(simbolos, puntos, tramos, ecd_x, ecd_y):
    
    x = sym.symbols('x')
    intervalos_x, intervalos_y = dividirEnIntervalos(ecd_x,ecd_y,puntos)

    condiciones = []
    for i in range(len(tramos)-1):
        condiciones.extend(get_condiciones_iniciales(tramos[i],tramos[i+1],puntos[i]))
    
    #print(condiciones)
    respuesta = sym.solvers.solve(condiciones, simbolos)
    #print(respuesta)
    

    for t in range (len(tramos)):
        for s in simbolos:
            try:
                tramos[t] = tramos[t].subs(s,respuesta[s])
            except KeyError:
                pass
            except AttributeError:
                pass
    
    coef = set()
    phi = []
    for t in tramos:
        if isinstance(t,int):
            phi.append(dict())
            continue

        phi2 = dict()
        for p in t.free_symbols:
            if p == x:
                continue
            coef.add(p)
            aux = t
            for p2 in t.free_symbols:
                if p2!=p and p2!=x:
                    aux = aux.subs(p2,0)
            phi2[p]=aux.subs(p,1)
    
        phi.append(phi2)
    
    #A,b = calcularSistemaMinimosCuadrados1(coef,phi,ecd_x,ecd_y,puntos)

    A,b = calcularSistemaMinimosCuadrados2(coef,phi,intervalos_x,intervalos_y)


    val_coef = np.linalg.solve(np.array(A),np.array(b))
    coef = list(coef)
    for t in range (len(tramos)):
        for c in range(len(coef)):
            try:
                tramos[t] = tramos[t].subs(coef[c],val_coef[c])
            except KeyError:
                pass
            except AttributeError:
                pass

    graficarFuncion(tramos,puntos,intervalos_x)
    return tramos
    

def calcularCalidadAjuste1(ecd_x,ecd_y,puntos,tramos):

    x = sym.symbols('x')
    t = 0
    sum = 0
    for i in range(len(ecd_x)):
        if t<len(puntos) and ecd_x[i] == puntos[t]:
            t+=1
        if isinstance(tramos[t],int):
            sum+=ecd_y[i]**2
        else:
            sum+= (ecd_y[i] - tramos[t].subs(x,ecd_x[i]))**2
      
    return sum**(1/2)

def calcularCalidadAjuste2(ecd_x,ecd_y,puntos,tramos):

    x = sym.symbols('x')
    t = 0
    val = []
    for i in range(len(ecd_x)):
        if t<len(puntos) and ecd_x[i] == puntos[t]:
            t+=1
        if isinstance(tramos[t],int):
            val.append(abs(ecd_y[i]))
        else:
            val.append(abs(ecd_y[i] - tramos[t].subs(x,ecd_x[i])))
      
    return max(val)
        

if __name__ == '__main__':
    '''prueba 1'''
    x,a2,b2,c2,a3,b3,c3,d3,a4,b4,c4,a5,b5,c5,d5,a6,b6,c6,a7,b7,c7,d7 = sym.symbols('x a2 b2 c2 a3 b3 c3 d3 a4 b4 c4 a5 b5 c5 d5 a6 b6 c6 a7 b7 c7 d7')
    puntos = [12.5, 24.5, 30, 39, 57.75, 66.5, 90]
    simbolos = [a2,b2,c2,a3,b3,c3,d3,a4,b4,c4,a5,b5,c5,d5,a6,b6,c6,a7,b7,c7,d7]
    tramos = [
    0,
    a2*x**2 + b2*x + c2,
    a3*x**3 + b3*x**2 + c3*x + d3,
    a4*x**2 + b4*x + c4,
    a5*x**3 + b5*x**2 + c5*x + d5,
    a6*x**2 + b6*x + c6,
    a7*x**3 + b7*x**2 + c7*x + d7,
    0
    ]
    #minimos_cuadrados(simbolos,puntos,tramos, ecd_x, ecd_y)

    '''prueba 2'''
    x,a2,b2,c2,a3,b3,c3,d3,a4,b4,c4,a5,b5,c5,d5,a6,b6,c6,a7,b7,c7,d7,d2,a8,b8,c8,d8,a9,b9,c9,d9 = sym.symbols('x a2 b2 c2 a3 b3 c3 d3 a4 b4 c4 a5 b5 c5 d5 a6 b6 c6 a7 b7 c7 d7 d2 a8 b8 c8 d8 a9 b9 c9 d9')
    puntos = [12.5, 24.5, 30, 33, 35, 39, 57.75, 66.5, 90]
    simbolos = [a2,b2,c2,a3,b3,c3,d3,a4,b4,c4,a5,b5,c5,d5,a6,b6,c6,a7,b7,c7,d7,d2,a8,b8,c8,d8,a9,b9,c9,d9]
    tramos = [
    0,
    a2*x**3 + b2*x**2 + c2*x + d2,
    a3*x**3 + b3*x**2 + c3*x + d3,
    a8*x**3 + b8*x**2 + c8*x + d8,
    a4*x**2 + b4*x + c4,
    a9*x**3 + b9*x**2 + c9*x + d9,
    a5*x**3 + b5*x**2 + c5*x + d5,
    a6*x**2 + b6*x + c6,
    a7*x**3 + b7*x**2 + c7*x + d7,
    0
    ]
    #minimos_cuadrados(simbolos,puntos,tramos, ecd_x, ecd_y)
    
    '''prueba 3'''
    x,a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,a4,b4,c4,a5,b5,c5,d5,a6,b6,c6,d6,a7,b7,c7,a8,b8,c8,d8,a9,b9,c9,d9 = sym.symbols('x a1 b1 c1 d1 a2 b2 c2 d2 a3 b3 c3 d3 a4 b4 c4 a5 b5 c5 d5 a6 b6 c6 d6 a7 b7 c7 a8 b8 c8 d8 a9 b9 c9 d9')
    puntos = [12.5, 21.5, 30, 33, 35, 39, 61, 71.5, 79, 90]
    simbolos = [a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,a4,b4,c4,a5,b5,c5,d5,a6,b6,c6,d6,a7,b7,c7,a8,b8,c8,d8,a9,b9,c9,d9]
    tramos = [
    0,
    a1*x**3 + b1*x**2 + c1*x + d1,
    a2*x**3 + b2*x**2 + c2*x + d2,
    a3*x**3 + b3*x**2 + c3*x + d3,
    a4*x**2 + b4*x + c4,
    a5*x**3 + b5*x**2 + c5*x + d5,
    a6*x**3 + b6*x**2 + c6*x + d6,
    a7*x**2 + b7*x + c7,
    a8*x**3 + b8*x**2 + c8*x + d8,
    a9*x**3 + b9*x**2 + c9*x + d9,
    0
    ]
    tramos = minimos_cuadrados(simbolos,puntos,tramos, ecd_x, ecd_y)

    # el tramo 1 es el que contiene al maximo de la P
    max_p = tramos[1].subs(x,sym.solvers.solve(sym.diff(tramos[1],x),x)[1])
    print('Maximo de la onda P =', max_p)
    print('Minimo de la onda P =',0)
    print('Amplitud de la onda P =',max_p)

    # el tramo 3 es el que contiene al mínimo de la QRS
    min_qrs = tramos[3].subs(x,sym.solvers.solve(sym.diff(tramos[3],x),x)[0])
    print('\nMinimo de la onda QRS =', min_qrs)

    # el tramo 4 es el que contiene al pico de la QRS
    max_qrs = tramos[4].subs(x,sym.solvers.solve(sym.diff(tramos[4],x),x)[0])
    print('Maximo de la onda QRS =', max_qrs)
    print('Amplitud de la onda QRS =',max_qrs-min_qrs)

   # el tramo 9 es el que contiene al mínimo de la T
    min_t = tramos[9].subs(x,sym.solvers.solve(sym.diff(tramos[9],x),x)[0])
    print('\nMinimo de la onda T =', min_t)

    # el tramo 7 es el que contiene al pico de la T
    max_t = tramos[7].subs(x,sym.solvers.solve(sym.diff(tramos[7],x),x)[0])
    print('Maximo de la onda T =', max_t)
    print('Amplitud de la onda T =',max_t-min_t)
    
    print('\nNorma 2 =',calcularCalidadAjuste1(ecd_x,ecd_y,puntos,tramos))

    print('\nNorma infinito = ',calcularCalidadAjuste2(ecd_x,ecd_y,puntos,tramos))
    
