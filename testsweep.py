import harmInstrI_SR860_update as harmInstrI
import matplotlib
matplotlib.rcParams['backend'] = 'wxagg'
import matplotlib.pyplot as plt
import numpy as np

lockin = harmInstrI.SR860()
lockin.filter = "10 ms"
X1, Y1, V1 = lockin.measureXYV(50, -50, 5)
X2, Y2, V2 = lockin.measureXYV(-50, 50, 5)

X = np.hstack((X1, X2))*1e12
Y = np.hstack((Y1, Y2))*1e12
V = np.hstack((V1, V2))

plt.plot(V, X)
plt.plot(V, Y)
plt.show()