import numpy as np
import json
import matplotlib

matplotlib.rcParams['backend'] = 'wxagg'

import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 1.5

fn = 'G:/My Drive/Data/boa/erogers2019/DOPS/20190723_006.dat'

d = json.load(open(fn, "r"))
Vac = 0.075*np.sqrt(2)
t = np.array([d['data'][j]['time'] for j in range(len(d['data']))])
capacitance = np.array([d['data'][j]['C'] for j in range(len(d['data']))])
h3x = np.array([d['data'][j]['harm3X'] for j in range(len(d['data']))])
h3y = np.array([d['data'][j]['harm3Y'] for j in range(len(d['data']))])
h3 = np.sqrt(h3x**2 + h3y**2)
v0 = np.array(d['pvals']['V0'])
v0err = np.array(d['perrs']['V0'])
pt = np.arange(len(d['data']))
slope = d['pvals']['b']
slope_err = d['perrs']['b']
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']

#plt.errorbar(pt, v0, v0err, fmt='-')
plt.plot(pt, h3/slope/Vac)
plt.xlabel('Time (s)')
plt.ylabel('Voltage Offset (mV)')
plt.show()