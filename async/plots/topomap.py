from asyncio.windows_events import NULL
import zlib
from eeg_positions import get_elec_coords, plot_coords
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CloughTocher2DInterpolator

#plot do escalpo e dos eletrodos em suas posições

RADIUS_INNER_CONTOUR = 0.72658518
RADIUS_OUTER_CONTOUR = 1.0

coords = get_elec_coords(
    elec_names=["AF3", "AF4", "T7", "T8", "Pz"],
    drop_landmarks=False,
    dim="2d",
)

fig, ax = plot_coords(
    coords, scatter_kwargs=dict(color="black"), text_kwargs=dict(fontsize=10)
)

#END plot do escalpo e dos eletrodos em suas posições


posX = [0, 0, 0, 0, 0]
posY = [0, 0, 0, 0, 0]

for idx, row in coords.iterrows():
    posX[idx] = row["x"]#array que guarda a coordenada x dos eletrodos
    posY[idx] = row["y"]#array que guara a coordenada y dos eletrodos

pointDensity = 500 #quantidade de pontos do linspace tanto no eixo x quanto no y
packetAmmount = 260
lastEpochDatapack = 260 #Último datapack considerado no epoch
firstEpochDatapack = 1 #Primeiro datapack considerado no epoch
epochLenght = (lastEpochDatapack - firstEpochDatapack + 1) / 8
drop  = False
dropList = NULL

if(firstEpochDatapack > 1):
    dropListBegin = np.array(np.arange(1, firstEpochDatapack, 1, dtype = int))
    dropList = np.append(dropList, dropListBegin)
    drop = True

if(lastEpochDatapack < packetAmmount):
    dropListEnd = np.array(np.arange(lastEpochDatapack + 1, packetAmmount, 1, dtype = int))
    dropList = np.append(dropList, dropListEnd)
    drop = True

if(drop):
    np.delete(dropList, 0)


Current_subject_file = 'Felipe-18-05-22.csv'
Current_subject = 'Felipe_01'

# Separei os dados do CSV de acordo com a faixa de frequência que eles representam, de forma com que é fácil trocar entre as faixas analisadas, bastando 
# descomentar o 'z' referente à faixa desejada ao passo que comenta o 'z' da faixa que estava sendo analisada.
data = pd.read_csv('../leituras/leituras_18-05-22/' + Current_subject_file)

if(drop):
    zAlpha = [data['af3a'].drop(dropList).mean(), data['af4a'].drop(dropList).mean(), data['pxa'].drop(dropList).mean(), data['t8a'].drop(dropList).mean(), data['pza'].drop(dropList).mean()]
    zLowBeta = [data['af3lb'].drop(dropList).mean(), data['af4lb'].drop(dropList).mean(), data['pxlb'].drop(dropList).mean(), data['t8lb'].drop(dropList).mean(), data['pzlb'].drop(dropList).mean()]
    zHighBeta = [data['af3hb'].drop(dropList).mean(), data['af4hb'].drop(dropList).mean(), data['pxhb'].drop(dropList).mean(), data['t8hb'].drop(dropList).mean(), data['pzhb'].drop(dropList).mean()]
    zGamma = [data['af3g'].drop(dropList).mean(), data['af4g'].drop(dropList).mean(), data['pxg'].drop(dropList).mean(), data['t8g'].drop(dropList).mean(), data['pzg'].drop(dropList).mean()]
    zTheta = [data['af3t'].drop(dropList).mean(), data['af4t'].drop(dropList).mean(), data['pxt'].drop(dropList).mean(), data['t8t'].drop(dropList).mean(), data['pzt'].drop(dropList).mean()]
else:
    zAlpha = [data['af3a'].mean(), data['af4a'].mean(), data['pxa'].mean(), data['t8a'].mean(), data['pza'].mean()]
    zLowBeta = [data['af3lb'].mean(), data['af4lb'].mean(), data['pxlb'].mean(), data['t8lb'].mean(), data['pzlb'].mean()]
    zHighBeta = [data['af3hb'].mean(), data['af4hb'].mean(), data['pxhb'].mean(), data['t8hb'].mean(), data['pzhb'].mean()]
    zGamma = [data['af3g'].mean(), data['af4g'].mean(), data['pxg'].mean(), data['t8g'].mean(), data['pzg'].mean()]
    zTheta = [data['af3t'].mean(), data['af4t'].mean(), data['pxt'].mean(), data['t8t'].mean(), data['pzt'].mean()]
    
#z = [zAlpha[0], zAlpha[1], zAlpha[2], zAlpha[3], zAlpha[4]]
#z = [zLowBeta[1], zLowBeta[2], zLowBeta[3], zLowBeta[4], zLowBeta[5]]
#z = [zHighBeta[1], zHighBeta[2], zHighBeta[3], zHighBeta[4], zHighBeta[5]]
#z = [zGamma[1], zGamma[2], zGamma[3], zGamma[4], zGamma[5]]
#z = [zTheta[1], zTheta[2], zTheta[3], zTheta[4], zTheta[5]]

# linspaces com as coordenadas dos pontos que preencherão o gráfico
xi = np.linspace(-1, 1, pointDensity)
yi = np.linspace(-1, 1, pointDensity)


# Para que a figura do escalpo seja completamente preenchida pelo heatmap,
# é necessário adicionar alguns outros pontos, em especial nas bordas da
# figura. Desta forma, foram adicionados 36 pontos igualmente espaçados 
# nas bordas do círculo do topomap.
newX = np.array(np.zeros((36, 1), float))
newY = np.array(np.zeros((36, 1), float))
newZa = np.array(np.zeros((36, 1), float))
newZlb = np.array(np.zeros((36, 1), float))
newZhb = np.array(np.zeros((36, 1), float))
newZg = np.array(np.zeros((36, 1), float))
newZt = np.array(np.zeros((36, 1), float))
distance = np.array(np.zeros((5, 1), float))
aux1 = 0
aux2 = 1

# Para calular o valor da potência nestes pontos adicionados, foi feita uma lógica de 
# interpolação baseada na distância entre os pontos adicionados e os eletrodos, sendo
# selecionados para a interpolação os dois eletrodos mais próximos.
for i in range(36):
    angle = (np.pi/18)*i
    newX[i] = np.cos(angle)
    newY[i] = np.sin(angle)
    distance[0]= np.sqrt((posX[0] - newX[i])**2 + (posY[0] - newY[i])**2)
    distance[1]= np.sqrt((posX[1] - newX[i])**2 + (posY[1] - newY[i])**2)
    distance[2] = np.sqrt((posX[2] - newX[i])**2 + (posY[2] - newY[i])**2)
    distance[3] = np.sqrt((posX[3] - newX[i])**2 + (posY[3] - newY[i])**2)
    distance[4] = np.sqrt((posX[4] - newX[i])**2 + (posY[4] - newY[i])**2)
    
    minDistance = distance[0]
    minDistance2 = distance[1]
    for j in range(3):
        if distance[j+2] < minDistance:
            minDistance, minDistance2 = distance[j+2], minDistance
            aux1, aux2 = j+2, aux1
        elif distance[j+2] < minDistance2:
            minDistance2 = distance[j+2]
            aux2 = j+2

    newZa[i] = ((1 - (distance[aux1] / 2)) * zAlpha[aux1])*0.75 + (((1 - (distance[aux2] / 2)) * zAlpha[aux2])*0.25)#lógica de interpolação
    newZlb[i] = ((1 - (distance[aux1] / 2)) * zLowBeta[aux1])*0.75 + (((1 - (distance[aux2] / 2)) * zLowBeta[aux2])*0.25)#lógica de interpolação
    newZhb[i] = ((1 - (distance[aux1] / 2)) * zHighBeta[aux1])*0.75 + (((1 - (distance[aux2] / 2)) * zHighBeta[aux2])*0.25)#lógica de interpolação
    newZg[i] = ((1 - (distance[aux1] / 2)) * zGamma[aux1])*0.75 + (((1 - (distance[aux2] / 2)) * zGamma[aux2])*0.25)#lógica de interpolação
    newZt[i] = ((1 - (distance[aux1] / 2)) * zTheta[aux1])*0.75 + (((1 - (distance[aux2] / 2)) * zTheta[aux2])*0.25)#lógica de interpolação

# FIM da adição de pontos para interpolação

#adição dos novos pontos gerados e interpolados nas estruturas que armazenam os dados que serão interpolados no topomap
posX  = np.append(posX, newX)
posY = np.append(posY, newY)
zAlpha = np.append(zAlpha, newZa)
zLowBeta = np.append(zLowBeta, newZlb)
zHighBeta = np.append(zHighBeta, newZhb)
zGamma = np.append(zGamma, newZg)
zTheta = np.append(zTheta, newZt)

X, Y = np.meshgrid(xi,yi) # "suaviza" as curvas no heatmap
interpAlpha = CloughTocher2DInterpolator(list(zip(posX, posY)), zAlpha) #interpolador das ondas alpha
Za = interpAlpha(X, Y)
interpLowBeta = CloughTocher2DInterpolator(list(zip(posX, posY)), zLowBeta) #interpolador das ondas low beta
Zlb = interpLowBeta(X, Y)
interpHighBeta = CloughTocher2DInterpolator(list(zip(posX, posY)), zHighBeta) #interpolador das ondas high beta
Zhb = interpHighBeta(X, Y)
interpGamma = CloughTocher2DInterpolator(list(zip(posX, posY)), zGamma) #interpolador das ondas gamma
Zg = interpGamma(X, Y)
interpTheta = CloughTocher2DInterpolator(list(zip(posX, posY)), zTheta) #interpolador das ondas theta
Zt = interpTheta(X, Y)


#lógica de plotagem
plt.pcolormesh(X, Y, Za, cmap = 'jet' , shading='auto')
plt.legend()
cbar = plt.colorbar()
cbar.set_label('uV²/Hz', rotation=90)

plt.axis("equal")
ax.axis("off")
plt.savefig(Current_subject + '_alpha.png')
cbar.remove()


plt.pcolormesh(X, Y, Zlb, cmap = 'jet' , shading='auto')
plt.legend()
cbar = plt.colorbar()
cbar.set_label('uV²/Hz', rotation=90)

plt.axis("equal")
ax.axis("off")
plt.savefig(Current_subject + '_lowBeta.png')
cbar.remove()

plt.pcolormesh(X, Y, Zhb, cmap = 'jet' , shading='auto')
plt.legend()
cbar = plt.colorbar()
cbar.set_label('uV²/Hz', rotation=90)

plt.axis("equal")
ax.axis("off")
plt.savefig(Current_subject + '_highBeta.png')
cbar.remove()

plt.pcolormesh(X, Y, Zg, cmap = 'jet' , shading='auto')
plt.legend()
cbar = plt.colorbar()
cbar.set_label('uV²/Hz', rotation=90)

plt.axis("equal")
ax.axis("off")
plt.savefig(Current_subject + '_gamma.png')
cbar.remove()

plt.pcolormesh(X, Y, Zg, cmap = 'jet' , shading='auto')
plt.legend()
cbar = plt.colorbar()
cbar.set_label('uV²/Hz', rotation=90)

plt.axis("equal")
ax.axis("off")
plt.savefig(Current_subject + '_theta.png')