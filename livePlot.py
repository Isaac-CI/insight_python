# -*- coding: utf-8 -*-
from itertools import count
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

plt.style.use('fivethirtyeight')

x_vals = []
y_vals = []

index = count()


def animate(i):
    data = pd.read_csv('data.csv')
    x = data['itr']
    y1 = data['af3a']
    y2 = data['af3hb']
    y3 = data['af3t']
    y4 = data['af4a']
    y5 = data['af4hb']
    y6 = data['af4t']
    #y7 = data['pxa']
    #y8 = data['pxhb']
    #y9 = data['pxt']
    
    plt.cla()

    for i in x:
        if i>600 and i<=1044: #Momento aprox do fim da calibragem
            if i % 37 == 0: #Tempo aprox do intervalo das imagens (n primo)
                if i % 2 == 0: #Definicao da cor
                    j = i
                    plt.axvspan(j, j+37, facecolor='b', alpha=0.3)
                else:
                    j = i
                    plt.axvspan(j, j+37, facecolor='r', alpha=0.3)
    #plt.axvspan(0, 672, facecolor='b', alpha=0.3)
    #plt.axvspan(672, i+1042, facecolor='g', alpha=0.3)

    plt.plot(x, y1, label='AF3Alpha')
    plt.plot(x, y2, label='AF3Hbeta')
    plt.plot(x, y3, label='AF3Theta')
    plt.plot(x, y4, label='Af4Alpha')
    plt.plot(x, y5, label='Af4HBeta')
    plt.plot(x, y6, label='Af4Theta')
    #plt.plot(x, y7, 'firebrick', label='PxAlpha')
    #plt.plot(x, y8, 'cyan', label='PxHBeta')
    #plt.plot(x, y9, 'black', label='PxTheta')

    plt.legend(loc='upper left')
    plt.tight_layout()
    #for i in x>672:



ani = FuncAnimation(plt.gcf(), animate, interval=1000)

plt.tight_layout()
plt.show()
