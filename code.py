import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.fft as sc

def hacerHorizontal(x,fx):
    recta_y = [] 
    ejex = []

    for i in range (1 , 409): 
        recta_y.append((i-1)*((fx[0]-fx[-1])/(x[0]-x[-1]))+fx[0])

    for i in range (0,fx.__len__()):
        fx[i] -= recta_y[i]
        ejex.append(0)

    return fxf

def primeraDerivada(f,paso):
    '''
    convolucion con d[n+1]-d[n-1] (aproximacion de la primera derivada)
    '''
    convolucion = np.convolve(f,[1,0,-1])
    for x in convolucion:
        x = x/(2*paso)

    return convolucion

def segundaDerivada(f,paso):
    '''
    convolucion con d[n+1]-2d[n]+d[n-1] (aproximacion de la segunda derivada)
    '''
    convolucion = np.convolve(f,[1,-2,1])
    for x in convolucion:
        x = x/(paso*paso)

    return convolucion

def hallarOndas(f, paso):
    derivada_1 = primeraDerivada(f, paso)
    criticos = []
    # Se consideran puntos críticos todos los puntos que sean 0, o aquellos que sean negativos y el siguiente a él positivo (y viceversa)
    for i in range (1, derivada_1.__len__()-1):
        if derivada_1[i] ==  0 or (derivada_1[i] < 0 and (derivada_1[i+1] > 0 or derivada_1[i-1] >0)) or (derivada_1[i] > 0 and (derivada_1[i+1] < 0 or derivada_1[i-1] < 0)): 
                criticos.append(i)
    
    # Los posibles puntos críticos que se encuentren sobre la onda QRS serán aquellos que la diferencia absoluta entre 2 sea mayor a 20 unidades.
    posible_qrs_y = []
    posible_qrs_x = []
    for x in criticos:
        for i in range (1,5):
            if x+i < derivada_1.__len__() and x-i > 0:
                if abs(derivada_1[x+i] - derivada_1[x]) > 20 and abs(derivada_1[x-i] - derivada_1[x]) > 20:
                    posible_qrs_x.append(x)
                    posible_qrs_y.append(derivada_1[x])

    qrs_x = []
    qrs_y = []
    # Finalmente, de los puntos críticos obtenidos solo se guarda aquellos (uno por periodo de la QRS) más cercanos a cero.
    for x in posible_qrs_x:
        posibles_x = []
        posibles_y = []
        for i in range (0, posible_qrs_x.__len__()):
            if abs(x - posible_qrs_x[i]) < 10:
                posibles_x.append(i)
                posibles_y.append(abs(posible_qrs_y[i]))
        indice = posible_qrs_y.index(min(posibles_y))
        try:
            qrs_y.index(posible_qrs_y[indice])
        except:
            qrs_x.append(posible_qrs_x[indice])
            qrs_y.append(posible_qrs_y[indice])
        
    solo_QRS = []
    solo_P = []
    solo_T = []  

    i = 0
    for x in qrs_x:
        while i<x-10:
            solo_QRS.append(0)
            solo_T.append(0)
            if i > x-40:
                solo_P.append(f[i])
            else:
                solo_P.append(0)
            i+=1
        for j in range(i,x+11):
            solo_QRS.append(f[j])
            solo_P.append(0)
            solo_T.append(0)
        i = x+11

        for j in range(i,i+60):
            solo_T.append(f[j])
            solo_QRS.append(0)
            solo_P.append(0)

        i=i+60
        
    return solo_QRS, solo_P, solo_T

def mostrarEspectroDeFrecuencias(f, duracion, titulo):
    n2 = f.__len__()
    t2 = np.linspace(0, duracion, n2) # Intervalo de tiempo en segundos
    dt2 = t2[1] - t2[0]
    y2 = f
    Y2 = sc.fft(f) 
    frq2 = sc.fftfreq(n2, dt2)
    fig = plt.figure(figsize=(6, 8))
    ax1 = fig.add_subplot(211)
    ax1.plot(t2, y2)
    plt.title('Espectro de Frecuencias ' + titulo)
    ax1.set_xlabel('Tiempo (s)')
    ax1.set_ylabel('$y_2(t)$')
    ax2 = fig.add_subplot(212)
    im = []
    for i in Y2:
        im.append(np.sqrt(i.real**2+i.imag**2))
    ax2.vlines(frq2, 0, im)
    plt.xlabel('Frecuencia (Hz)')
    plt.ylabel('Mod($Y_2$)')
    plt.show()
    

def filtrarOnda(f, porcentaje):
    F = sc.fft(f)
    im = []
    for i in F:
        im.append(np.sqrt(i.real**2+i.imag**2))
    maximo_modulo = max(im)
    for i in range (0, im.__len__()):
        if(im[i] < porcentaje*maximo_modulo): F[i] = 0
    return sc.ifft(F)

if __name__ == "__main__":

    #Tomar datos del archivo cardio.csv
    data = pd.read_csv("cardio.csv",",",header=None)
    ecd_x = data[0].tolist()
    ecd_y = data[1].tolist()

    #Modificar horizontalidad de la funcion
    ecd_y = hacerHorizontal(ecd_x,ecd_y)
    #Graficar
    plt.title("Senal Horizontal")
    plt.plot(ecd_x, ecd_y)
    plt.show()

    #Hallar Derivadas
    derivada_1 = primeraDerivada(ecd_y,180/408)
    derivada_2 = segundaDerivada(ecd_y,180/408)
    #Graficar
    plt.title("Primer Derivada")
    plt.plot(np.arange(0,derivada_1.__len__(),1),derivada_1)
    plt.show()
    plt.title("Segunda Derivada")
    plt.plot(np.arange(0,derivada_2.__len__(),1),derivada_2)
    plt.show()

    #Separar Ondas P, QRS y T
    solo_QRS, solo_P, solo_T = hallarOndas(ecd_y, 180/408)
    plt.title("Ondas Separadas")
    plt.plot(np.arange(0,solo_QRS.__len__(),1),solo_QRS, label="Onda QRS")
    plt.plot(np.arange(0,solo_P.__len__(),1),solo_P, label="Onda P")
    plt.plot(np.arange(0,solo_T.__len__(),1),solo_T, label="Onda T")
    plt.legend(loc="best")
    plt.show()

    #Graficar Espectros de Frecuencia Ondas
    mostrarEspectroDeFrecuencias(solo_QRS,180, "Onda QRS")
    mostrarEspectroDeFrecuencias(solo_P,180, "Onda P")
    mostrarEspectroDeFrecuencias(solo_T,180, "Onda T")

    #Filtrar Ondas
    QRS_filtrada = filtrarOnda(solo_QRS, 0.12)
    P_filtrada = filtrarOnda(solo_P, 0.4)
    T_filtrada = filtrarOnda(solo_T, 0.12)

    #Reconstruir Onda
    reconstruccion = QRS_filtrada + P_filtrada + T_filtrada
    #Graficar
    plt.title("Original VS Reconstruida Filtrada")
    plt.plot(ecd_x, ecd_y,label="Original")
    plt.plot(np.arange(0, reconstruccion.__len__(), 1), reconstruccion, label="Reconstruida Filtrada")
    plt.legend(loc="best")
    plt.show()

    #Filtrar Senal Original
    ecd_y_filtrada = filtrarOnda(ecd_y,0.05)
    #Comparaciones Graficas
    plt.title("Original VS Original Filtrada")
    plt.plot(ecd_x,ecd_y, label="Original")
    plt.plot(ecd_x, ecd_y_filtrada, label = "Original Filtrada")
    plt.legend(loc="best")
    plt.show()

    plt.title("Original Filtrada VS Reconstruida Filtrada")
    plt.plot(ecd_x, ecd_y_filtrada, label = "Original Filtrada")
    plt.plot(np.arange(0, reconstruccion.__len__(), 1), reconstruccion, label="Reconstruida Filtrada")
    plt.legend(loc="best")
    plt.show()

    #Graficar Espectro de Frecuencia
    mostrarEspectroDeFrecuencias(reconstruccion, 180, "Senal Reconstruida Filtrada")
    mostrarEspectroDeFrecuencias(ecd_y_filtrada,180, "Senal Original Filtrada")

    plt.ion()

