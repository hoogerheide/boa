import numpy as np
import json
import matplotlib

matplotlib.rcParams['backend'] = 'wxagg'

import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 1.5

fn = 'DOPCDOPS/20190723_005.dat'

d = json.load(open(fn, "r"))

t = np.array([d['data'][j]['time'] for j in range(len(d['data']))])
capacitance = np.array([d['data'][j]['C'] for j in range(len(d['data']))])
v0 = np.array(d['pvals']['V0'])
v0err = np.array(d['perrs']['V0'])
pt = np.arange(len(d['data']))
slope = d['pvals']['b']
slope_err = d['perrs']['b']
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']

plt.errorbar(pt, v0, v0err, fmt='-')
plt.xlabel('Time (s)')
plt.ylabel('Voltage Offset (mV)')
plt.show()