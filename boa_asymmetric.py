import numpy as np
import json
import matplotlib

matplotlib.rcParams['backend'] = 'wxagg'

import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 1.5

#trial 1
fn1 = 'DOPCDOPS/20190723_003.dat'
d1 = json.load(open(fn1, "r"))
t1 = np.array([d1['data'][j]['time'] for j in range(len(d1['data']))])
pt1 = np.arange(len(d1['data']))
v01 = np.array(d1['pvals']['V0'])
v0err1 = np.array(d1['perrs']['V0'])
slope1 = d1['pvals']['b']
slope_err1 = d1['perrs']['b']
offset1 = np.mean(v01[214:253])
stdev1 = np.std(v01[214:253])
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']

#trial 2
fn2 = 'DOPCDOPS/20190723_004.dat'
d2 = json.load(open(fn2, "r"))
t2 = np.array([d2['data'][j]['time'] for j in range(len(d2['data']))])
pt2 = np.arange(len(d2['data']))
v02 = np.array(d2['pvals']['V0'])
v0err2 = np.array(d2['perrs']['V0'])
slope2 = d2['pvals']['b']
slope_err2 = d2['perrs']['b']
offset2 = np.mean(v02[125:150])
stdev2 = np.std(v02[125:150])
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']

#trial 3
fn3 = 'DOPCDOPS/20190723_005.dat'
d3 = json.load(open(fn3, "r"))
t3 = np.array([d3['data'][j]['time'] for j in range(len(d3['data']))])
pt3 = np.arange(len(d3['data']))
v03 = np.array(d3['pvals']['V0'])
v0err3 = np.array(d3['perrs']['V0'])
slope3 = d3['pvals']['b']
slope_err3 = d3['perrs']['b']
offset3 = np.mean(v03[125:155])
stdev3 = np.std(v03[125:155])
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']
#PLOT THIS

stdev = np.std([stdev1, stdev2, stdev3])
offset = np.mean([offset1, offset2, offset3])
print(stdev, offset, stdev1, stdev2, stdev3, offset1, offset2, offset3)