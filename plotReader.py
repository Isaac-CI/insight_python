import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

filePath = 'async\leituras\Leituras_18-05-22\Jhonathan-18-05-22.csv'

plt.style.use('fivethirtyeight')

x = []
y = []
  
data = pd.read_csv(filePath)
x = data['itr']
lastPacket = data['itr'].iloc[-1]
y1 = data['af3a'] / data['af3a'].mean()
y2 = data['af3hb'] / data['af3hb'].mean()
y3 = data['af3t'] / data['af3t'].mean()
y4 = data['af4a'] / data['af4a'].mean()
y5 = data['af4hb'] / data['af4hb'].mean()
y6 = data['af4t'] / data['af4t'].mean()

ordenadas = np.zeros((lastPacket,), dtype=float)

for i in data['itr']:
    ordenadas[i - 1] = i / 8

  
plt.tight_layout()
plt.plot(ordenadas, y1, label='AF3Alpha')
plt.plot(ordenadas, y2, label='AF3Hbeta')
#plt.plot(ordenadas, y3, label='AF3Theta')
plt.plot(ordenadas, y4, label='Af4Alpha')
plt.plot(ordenadas, y5, label='Af4HBeta')
#plt.plot(ordenadas, y6, label='Af4Theta')
plt.legend(loc='upper left')
plt.plot()
plt.show()