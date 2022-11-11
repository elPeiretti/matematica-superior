import numpy as np
import sympy as sym

def derivadasParciales(f):
        return (sym.diff(f,x), sym.diff(f,y), sym.diff(f,z))

def evaluar(f,P):
    return f.subs([(x,P[0]),(y,P[1]),(z,P[2])])

def newton_raphson(x0, f, max_iter, tol):
    iter = 0
    ans = []
    xi = x0
    t = sym.symbols('t')
    dfdt = sym.diff(f,t)

    f_evaluada = f.subs(t,x0)
    
    ans.append(f_evaluada)

    # codigo = 0 -> retorna por condicion de tolerancia
    # codigo = -1 -> retorna por divergencia del metodo
    # codigo = 1 -> retorna por condicion de maximas iteraciones

    while iter<=max_iter and abs(ans[-1])>tol:
        
        # intentamos calcular el xi haciendo uso de fracciones, pero en cada iteracion eran cada vez más grandes
        # y eso hacía que el método se vuelva muy lento
        xi = float(xi - (f_evaluada/dfdt.subs(t,xi)))
        f_evaluada = f.subs(t,xi)
        ans.append(f_evaluada)
        
        if iter>1 and abs(ans[-2]) < abs(ans[-1]): #diverge el metodo
            return xi, ans[-1], -1
        iter += 1

    if iter > max_iter:
        return xi, ans[-1], 1
    
    return xi, ans[-1], 0

    

# f y g son las funciones cuyas superficies de nivel seran paralelas
# a la recta r que pasa por el punto P, y h es la funcion cuya superficie
# de nivel se utiliza para intersectar con la recta r y obtener el punto Q
def obtenerPuntoInterseccion(f, g, h, P, tolerancia, max_iter):

    fx,fy,fz = derivadasParciales(f)
    gx,gy,gz = derivadasParciales(g)

    derivf_evaluada = [float(evaluar(fx,P)), float(evaluar(fy,P)), float(evaluar(fz,P))]
    derivg_evaluada = [float(evaluar(gx,P)), float(evaluar(gy,P)), float(evaluar(gz,P))]

    director1 = np.cross(derivf_evaluada,derivg_evaluada)
    t = sym.symbols('t')
    r = (director1[0]*t + P[0], director1[1]*t + P[1], director1[2]*t + P[2])
    #print(r)
    intersec_h_r1 = evaluar(h,r)
    #print(intersec_h_r1)

    # Obtención del punto Q
    # ==========================================================================================================
    # Pensamos en observar graficamente los puntos de las intersecciones para calcular el xn inicial para newton
    # pero nos pareció muy poco "programatico".
    # Pensamos tambien en aplicar biseccion o regula falsi para obtener un xn que estamos seguro que va a converger, pero no tenemos forma de determinar el intervalo [izq,der] para dichos metodos.
    # Por lo tanto, se decidió en tomar un punto inicial arbitrario para comenzar el método y utilizar este para todos los casos
    # ==========================================================================================================
    
    # los argumentos son: t0, f, max_iteraciones, tolerancia
    tn, yn, codigo = newton_raphson(1,intersec_h_r1,max_iter,tolerancia*0.01)
    if codigo != 0:
        return codigo

    Q = (r[0].subs(t,tn), r[1].subs(t,tn), r[2].subs(t,tn))
    return Q
    
def validarPunto(P):

    if not(type(P) is int):
        return "ok"

    if P == -1:
        return "EL METODO N-R DIVERGE"
    if P == 1:
        return "N-R FINALIZO POR CANTIDAD DE ITERACIONES"

def reboteMS(f1, f2, f3, Pn, tolerancia, max_iter):
# max_iter será la cantidad de iteraciones máximas para el método N-R

    Q1 = obtenerPuntoInterseccion(f1,f2,f3,Pn,tolerancia,max_iter)
    msg = validarPunto(Q1)
    if msg != "ok":
        return msg

    Q2 = obtenerPuntoInterseccion(f2,f3,f1,Q1,tolerancia,max_iter)
    msg = validarPunto(Q2)
    if msg != "ok":
        return msg

    Pn1 = obtenerPuntoInterseccion(f1,f3,f2,Q2,tolerancia,max_iter)
    msg = validarPunto(Pn1)
    if msg != "ok":
        return msg
    
    return Pn1

def dist(P1,P2):
    return (P1[0]-P2[0])**2 + (P1[1]-P2[1])**2 + (P1[2]-P2[2])**2

def reboteMS_iterativo(f1, f2, f3, Pn, tolerancia, max_iter):
    
    i=0
    puntos = []
    ans = [max(abs(float(evaluar(f1,Pn))), abs(float(evaluar(f2,Pn))), abs(float(evaluar(f3,Pn))))]

    while not(isinstance(Pn,str)) and i<=max_iter and ans[-1] > tolerancia:
        puntos.append(Pn)
        
        if i>1 and abs(ans[-2]) < abs(ans[-1]):
            print("El metodo reboteMS diverge")
            return

        Pn = reboteMS(f1,f2,f3, Pn, tolerancia, max_iter*2)
        if not isinstance(Pn,str):
            ans.append(max(abs(float(evaluar(f1,Pn))), abs(float(evaluar(f2,Pn))), abs(float(evaluar(f3,Pn)))))
        i+=1

    if isinstance(Pn,str):
        print("El método falló en la iteración", i)
        print("La causa del fallo es:", Pn)

    print("iteraciones:",i)
    return Pn

def newton_raphson_sistema(f, g, h, Pn, max_iter, tolerancia):
    jacobiano = []
    jacobiano.append(np.asarray(derivadasParciales(f)))
    jacobiano.append(np.asarray(derivadasParciales(g)))
    jacobiano.append(np.asarray(derivadasParciales(h)))
    iter = 0
    ans = [max(float(abs(evaluar(f,Pn))), abs(float(evaluar(g,Pn))), abs(float(evaluar(h,Pn))))]

    # codigo = 0 -> retorna por condicion de tolerancia
    # codigo = 1 -> retorna por condicion de maximas iteraciones
    # codigo = -1 -> retorna por condicion de divergencia

    while max_iter>iter and ans[-1] > tolerancia:
        
        #evaluacion de las derivadas parciales en el punto
        jacobiano_evaluado = []
        for deriv in jacobiano:
            jacobiano_evaluado.append([float(evaluar(deriv[0],Pn)),float(evaluar(deriv[1],Pn)),float(evaluar(deriv[2],Pn))])
        #evaluacion de las funciones en el punto
        funciones_evaluadas = [[-1*float(evaluar(f,Pn))],[-1*float(evaluar(g,Pn))],[-1*float(evaluar(h,Pn))]]
        
        if np.linalg.det(jacobiano_evaluado) == 0:
            print("iteraciones:",iter)
            print("NO ES INVERTIBLE LA MATRIZ, NO SE PUEDE CONTINUAR")
            return

        inversa_jacobiano = np.matrix(jacobiano_evaluado)**-1
        
        delta =  inversa_jacobiano * funciones_evaluadas

        Pn = (Pn[0] + delta[0,0], Pn[1] + delta[1,0], Pn[2] + delta[2,0])
       
        ans.append(max(abs(float(evaluar(f,Pn))), float(abs(evaluar(g,Pn))), float(abs(evaluar(h,Pn)))))
        
        if iter > 0 and ans[-1] > ans[-2]:
            print("iteraciones:",iter)
            return Pn, ans, -1
        
        iter+=1
    print("iteraciones:",iter)

    
    if iter < max_iter:
        return Pn, ans, 0

    return Pn, ans, 1

if __name__  == "__main__":
    x,y,z = sym.symbols('x y z')
    f1 = sym.cos(x) + sym.cos(y) - z + 1
    f2 = x**2  + 3 * (y**2) - z
    f3 = x + y - z

    tol = 0.00000001
    it = 100
    P = (2,-0.5,1)

    import time
    print("=======REBOTE=====")
    inicio = time.time()
    Pn1 = reboteMS_iterativo(f1,f2,f3, P, tol, it)
    print(Pn1)
    print("tiempo [s]:",time.time()-inicio)

    print("=======NEWTON=====")
    inicio = time.time()
    Pn1, ans, codigo = newton_raphson_sistema(f1,f2,f3,P,it,tol)
    print("tiempo [s]:",time.time()-inicio)
    if(codigo == 0):
        print("Solucion al sistema:",Pn1)
    elif codigo == -1:
        print("El método diverge.",Pn1)
    else:
        print("El método finalizó por cantidad de iteraciones. Pn=",Pn1)
