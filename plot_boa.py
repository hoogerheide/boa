import numpy as np
import matplotlib

matplotlib.rcParams['backend'] = 'wxagg'

import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 1.5

fn = '20190717_001.txt'

d = np.loadtxt(fn)

plt.plot(d[:,0], d[:,2])
plt.show()