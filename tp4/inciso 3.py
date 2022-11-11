import numpy as np
import matplotlib.pyplot as plt
import time

#aproxima el siguiente paso
def runge_kutta4(m, paso):

    v = m[-1][1]
    w = m[-1][2]

    k1v = 3*(v+ w - (v**3)/3)
    k1w = -1/3*(w - 0.7 + 0.8*v)

    k2v = 3*( (v+k1v*paso/2) + (w+k1w*paso/2) - ((v+k1v*paso/2)**3)/3)
    k2w = -1/3*((w+k1w*paso/2) - 0.7 + 0.8*(v+k1v*paso/2))

    k3v = 3*( (v+k2v*paso/2) + (w+k2w*paso/2) - ((v+k2v*paso/2)**3)/3)
    k3w = -1/3*((w+k2w*paso/2) - 0.7 + 0.8*(v+k2v*paso/2))

    k4v = 3*( (v+k3v*paso) + (w+k3w*paso) - ((v+k3v*paso)**3)/3)
    k4w = -1/3*((w+k3w*paso) - 0.7 + 0.8*(v+k3v*paso))

    v_sig = v + 1/6*paso*(k1v + 2*k2v + 2*k3v + k4v)
    w_sig = w + 1/6*paso*(k1w + 2*k2w + 2*k3w + k4w)

    return np.append(m,[[m[-1][0]+paso, v_sig, w_sig, 3*(v_sig+ w_sig - (v_sig**3)/3), -1/3*(w_sig - 0.7 + 0.8*v_sig)]],axis=0)

''' [ 
      [t0   v0   w0  v'0  w'0]
      [t1   v1   w1  v'1  w'1]
      [t2   v2   w2  v'2  w'2]
      [t3   v3   w3  v'3  w'3]
      ......................  
    ]
'''

def metodo_pcs(m,paso):
    
    t = m[-1][0]

    #predictor
    vp = m[-3][1] + paso * (0.375*m[-3][3] + 1.125*m[-1][3])
    wp = m[-3][2] + paso * (0.375*m[-3][4] + 1.125*m[-1][4])
    dvp = 3*(vp+ wp - (vp**3)/3)
    dwp = -1/3*(wp - 0.7 + 0.8*vp)

    #corrector
    vc = m[-2][1] + 0.5*(paso/3)*(m[-2][3] + 4*m[-1][3] + dvp)
    wc = m[-2][2] + 0.5*(paso/3)*(m[-2][4] + 4*m[-1][4] + dwp)
    dvc = 3*(vc+ wc - (vc**3)/3)
    dwc = -1/3*(wc - 0.7 + 0.8*vc)

    m = np.append(m, [[m[-1][0]+paso*0.5, vc, wc, dvc, dwc]], axis=0)

    #supracorrector
    vsp = m[-4][1] + paso*( (11/12)*(m[-3][3]+(1/22)*m[-5][3]) + (1/24)*dvc)
    wsp = m[-4][2] + paso*( (11/12)*(m[-3][4]+(1/22)*m[-5][4]) + (1/24)*dwc)
    dvsp = 3*(vsp+ wsp - (vsp**3)/3)
    dwsp = -1/3*(wsp - 0.7 + 0.8*vsp)

    m[-2][1] = vsp
    m[-2][2] = wsp
    m[-2][3] = dvsp
    m[-2][4] = dwsp

    return m, (17/15)*(vc-vsp), (17/15)*(wc-wsp)

paso = 0.002
matriz = np.array([[0,0,0,0,7/30]])
v = []
w = []
t = []
diferencia_v = []
diferencia_w = []
richardson_v = []
richardson_w = []


start = time.time()

# 3 iteraciones de RK
for i in range(3):
    v.append(matriz[-1][1])
    w.append(matriz[-1][2])
    t.append(matriz[-1][0])
    matriz = runge_kutta4(matriz,0.5*paso)

# metodo PCS
for i in range(10000):
    v.append(matriz[-1][1])
    w.append(matriz[-1][2])
    t.append(matriz[-1][0])
    matriz, rich_v, rich_w = metodo_pcs(matriz,paso)
    richardson_v.append(rich_v)
    richardson_w.append(rich_w)


print("METODO PCS TARDO --- %s seconds ---" % (time.time() - start))

v1 = []
w1 = []
t1 = []
matriz = np.array([[0,0,0,0,7/30]])

start = time.time()

for i in range(10000):
    v1.append(matriz[-1][1])
    w1.append(matriz[-1][2])
    t1.append(matriz[-1][0])
    matriz = runge_kutta4(matriz,0.5*paso)

print("METODO R-K 4 TARDO --- %s seconds ---" % (time.time() - start))

# V y W en funcion del tiempo
plt.plot(t,v,label="v PCS")
plt.plot(t,w,label="w PCS")
plt.plot(t1,v1,label="v RK")
plt.plot(t1,w1,label="w RK")
plt.legend(loc='best')
plt.xlabel("Tiempo")
plt.ylabel("V | W")
plt.grid(True)
plt.show()

#estimador de richardson en funcion del tiempo
plt.scatter(t[3:],richardson_v,label="error estimado V")
plt.scatter(t[3:],richardson_w,label="error estimado W")
plt.xlabel("Tiempo")
plt.ylabel("Error estimado")
plt.legend(loc="best")
plt.grid(True)
plt.show()

# diagrama V-W
plt.plot(v,w,label="metodo PCS")
plt.plot(v1,w1,label="metodo R-K 4")
plt.xlabel("V")
plt.ylabel("W")
plt.legend(loc="best")
plt.grid(True)
plt.show()
