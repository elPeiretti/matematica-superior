from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

celeste = (154, 194, 231, 255)
rojo = (255,0,0,255)
naranja = (255,128,0,255)

def obtenerMalla(paso):
    # Paso*Paso = tamano del nodo

    im = Image.open("corazon.png")
    pix = im.load()

    h=int(np.ceil(im.size[0]/paso))
    w=int(np.ceil(im.size[1]/paso))

    nueva_imagen = np.zeros((h, w, 4), dtype=np.uint8)
    ans = []

    for i in range(0,im.size[0]-paso,paso):
        aux = []
        for j in range(0,im.size[1]-paso,paso):
            c=0
            r=0
            n=0
            for a in range(i,i+paso):
                for b in range(j,j+paso):
                    if pix[a,b] == celeste:
                        c+=1
                    elif pix[a,b] == naranja:
                        n+=1
                    else:
                        r+=1
            x = max(c,r,n)
            if x==c:
                nueva_imagen[int(i/paso)][int(j/paso)] = np.array(celeste)
                aux.append("celeste")
            elif x==n:
                nueva_imagen[int(i/paso)][int(j/paso)] = np.array(naranja)
                aux.append("naranja")
            else:
                nueva_imagen[int(i/paso)][int(j/paso)] = np.array(rojo)
                aux.append("rojo")
        ans.append(aux)
            

    #im2 = Image.fromarray(nueva_imagen,"RGBA")
    #im2.show()
    return ans

def obtenerMatrizTemperaturasIniciales(malla):
    ans = []
    for i in range(len(malla)):
        aux = []
        for j in range(len(malla[0])):
            if malla[i][j] == "celeste":
                        aux.append(-150)
            else:
                aux.append(36)
        ans.append(aux)
    
    return ans
    
def calcularTemperaturaSiguienteTiempo(malla, temperaturas, pasoEspacial, pasoTemporal):
    ans = []
    enfriado = True
    for i in range(len(malla)):
        aux = []
        for j in range(len(malla[0])):
            
            coef_k = 1
            if malla[i][j] == "naranja":
                coef_k = 0.8
            elif malla[i][j] == "rojo":
                coef_k=0.5

            temp = temperaturas[i][j]
            coef = coef_k*pasoTemporal/(pasoEspacial**2)

            #centro
            temp += coef*-4*temperaturas[i][j]

            #arriba
            if i==0:
                temp += coef*-150
            else:
                temp += coef*temperaturas[i-1][j]

            #abajo
            if i==len(malla)-1:
                temp += coef*-150
            else:
                temp += coef*temperaturas[i+1][j]

            #izq
            if j==0:
                temp += coef*-150
            else:
                temp += coef*temperaturas[i][j-1]

            #der
            if j==len(malla[0])-1:
                temp += coef*-150
            else:
                temp += coef*temperaturas[i][j+1]
            
            if temp >= 0:
                enfriado = False

            aux.append(temp)
    
        ans.append(aux)

    return ans, enfriado


pasoEspacial = 10
pasoTemporal = 0.0000001
conversor_px_metros = 0.00017

import time
start = time.time()

malla = obtenerMalla(pasoEspacial)
temperaturas = obtenerMatrizTemperaturasIniciales(malla)

fig,ax = plt.subplots(1,1)
imag = ax.imshow(temperaturas,cmap=plt.get_cmap('jet'),vmin=-200,vmax=100)
plt.colorbar(imag,ax=ax,label="Temperatura [°C]")

tiempo = 0
bajo_cero = False
while not bajo_cero:
    temperaturas, bajo_cero = calcularTemperaturaSiguienteTiempo(malla,temperaturas,pasoEspacial*conversor_px_metros,pasoTemporal)
    tiempo+=pasoTemporal
    print(tiempo)
    imag.set_data(temperaturas)
    fig.canvas.draw_idle()
    # descomentar la línea de abajo para visualizar el método
    #plt.pause(0.0001)

print("--- %s seconds ---" % (time.time() - start))
plt.show()
# h1 = paso temporal
# h = pixeles * metros/pixeles