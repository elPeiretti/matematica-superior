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

print("Ejemplo 0:")
metodo_newton_final(    
                        [1,2,3,4,5]         ,
                        [2,4,8,16,32]       , 
                        [[1,3,5],[3,5,-20]] ,
                        [[1,3],[-1,70]]     ,
                        [[1,3],[50,2]]             
                    )


print("\nEjemplo 1:")
#Ejemplo1R
metodo_newton_final(    
                        [-5,-2,3,8,16]            ,
                        [0,7,-5,-7,-9]            , 
                        [  [-5,3]  , [2,-2]   ]    ,
                        [  [3]  , [-5]   ]        ,
                        [  [3]  , [10]   ]             
                    )


print("\nEjemplo 2:")
#Ejemplo2R
metodo_newton_final(    
                        [5,10,30,20]        ,
                        [15,10,15,5]         , 
                        [  [5,10]  , [-2,-2]   ]    ,
                        [  [5]  , [-1]   ]        ,
                        [  []  , []   ]             
                    )

print("\nEjemplo 3:")
#Ejemplo3R
metodo_newton_final(    
                        [-5,-4,1,5]            ,
                        [5,2,1,-3]            , 
                        [  [-5,5,1]  , [-4,-3,1]   ]    ,
                        [  [1,5]  , [1,50]   ]        ,
                        [  [1]  , [-5]   ]             
                    )

print("\nEjemplo 4:")
#Ejemplo4
metodo_newton_final(
                    [-3,-2,-1,0,1],
                    [3,-2,1,0,-1],
                    [ [-3,-2,-1] , [1,1,1] ],
                    [ [-3,-2,-1] , [1,1,1] ],
                    [ [-3,-2,-1] , [1,1,1] ]
)


''' 
Los argumentos son:
*El vector de coordenadas x de los puntos por donde tiene que pasar el polinomio [x0,x1,x2,..,xn]
*El vector de coordenadas y de los puntos por donde tiene que pasar el polinomio [y0,y1,y2,...yn]
*Un vector con dos vectores:
    -El primero contiene las coordenadas x de la primera derivada 
    -El segundo los valores de las derivadas en las x dadas 
    [[x1,x3,x8],[13,4,1]] significa que f'(x1)=13, f'(x3)=4 y f'(x8)=1
*Un vector con dos vectores:
    -El primero contiene las coordenadas x de la segunda derivada 
    -El segundo los valores de las derivadas en las x dadas 
*Un vector con dos vectores:
    -El primero contiene las coordenadas x de la tercera derivada 
    -El segundo los valores de las derivadas en las x dadas 
*Si no se quiere condicionar alguna derivada, directamente hay que pasarle los vectores vacios [[],[]]
'''
