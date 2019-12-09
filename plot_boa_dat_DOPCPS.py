import numpy as np
import json
import matplotlib

matplotlib.rcParams['backend'] = 'wxagg'

import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 1.5

#trial 1
fn1 = 'DOPC/20190725_003.dat'
d1 = json.load(open(fn1, "r"))
t1 = np.array([d1['data'][j]['time'] for j in range(len(d1['data']))])
pt1 = np.arange(len(d1['data']))
v01 = np.array(d1['pvals']['V0'])
v0err1 = np.array(d1['perrs']['V0'])
slope1 = d1['pvals']['b']
slope_err1 = d1['perrs']['b']
offset1 = np.mean(v01[185:192])
means1 = np.mean(v01[200:206])-offset1, np.mean(v01[211:217])-offset1, np.mean(v01[227:232])-offset1, np.mean(v01[242:249])-offset1, np.mean(v01[258:265])-offset1, np.mean(v01[274:279])-offset1, np.mean(v01[283:289])-offset1
#means = np.array(np.mean(ploc[])
stdevs1 = np.std(v01[200:206]), np.std(v01[211:217]), np.std(v01[227:232]), np.std(v01[242:249]), np.std(v01[258:265]), np.std(v01[274:279]), np.std(v01[283:289])
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']
conc1 = np.array([0.07424, 0.3669, 0.7261, 3.559, 7.038, 34.49, 68.23])
# fig, ax = plt.subplots()
# x_vals = np.array([100, 500, 1000, 5000, 10000, 50000, 100000])

#trial 2
fn2 = 'DOPC/20190729_006.dat'
d2 = json.load(open(fn2, "r"))
t2 = np.array([d2['data'][j]['time'] for j in range(len(d2['data']))])
pt2 = np.arange(len(d2['data']))
v02 = np.array(d2['pvals']['V0'])
v0err2 = np.array(d2['perrs']['V0'])
slope2 = d2['pvals']['b']
slope_err2 = d2['perrs']['b']
offset2 = np.mean(v02[127:152])
means2 = np.mean(v02[162:171])-offset2, np.mean(v02[179:189])-offset2, np.mean(v02[211:218])-offset2, np.mean(v02[228:232])-offset2, np.mean(v02[240:248])-offset2, np.mean(v02[271:277])-offset2
#means = np.array(np.mean(ploc[])
stdevs2 = np.std(v02[162:171]), np.std(v02[179:189]), np.std(v02[211:218]), np.std(v02[228:232]), np.std(v02[240:248]), np.std(v02[271:277])
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']
conc2 = np.array([0.07366, 0.3640, 0.7205, 3.532, 6.986, 34.24])
# fig, ax = plt.subplots()
# x_vals = np.array([100, 500, 1000, 5000, 10000, 50000, 100000])

#trial 3
fn3 = 'DOPC/20190729_009.dat'
d3 = json.load(open(fn3, "r"))
t3 = np.array([d3['data'][j]['time'] for j in range(len(d3['data']))])
pt3 = np.arange(len(d3['data']))
v03 = np.array(d3['pvals']['V0'])
v0err3 = np.array(d3['perrs']['V0'])
slope3 = d3['pvals']['b']
slope_err3 = d3['perrs']['b']
offset3 = np.mean(v03[120:130])
means3 = np.mean(v03[139:144])-offset3, np.mean(v03[148:155])-offset3, np.mean(v03[161:166])-offset3, np.mean(v03[172:176])-offset3, np.mean(v03[222:228])-offset3, np.mean(v03[235:245])-offset3, np.mean(v03[271:282])-offset3
#means = np.array(np.mean(ploc[])
stdevs3 = np.std(v03[139:144]), np.std(v03[148:155]), np.std(v03[161:166]), np.std(v03[172:176]), np.std(v03[222:228]), np.std(v03[235:245]), np.std(v03[271:282])
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']
conc3 = np.array([0.07786, 0.3845, 0.7606, 3.725, 7.363, 36.05, 71.28])
# fig, ax = plt.subplots()
# x_vals = np.array([100, 500, 1000, 5000, 10000, 50000, 100000])
#graphs
plt.errorbar(conc1, means1, stdevs1, fmt='o')
plt.errorbar(conc2, means2, stdevs2, fmt='o')
plt.errorbar(conc3, means3, stdevs3, fmt='o')

#trial 1 PS
fn1 = 'DOPS/20190724_002.dat'
d1 = json.load(open(fn1, "r"))
t1 = np.array([d1['data'][j]['time'] for j in range(len(d1['data']))])
pt1 = np.arange(len(d1['data']))
v01 = np.array(d1['pvals']['V0'])
v0err1 = np.array(d1['perrs']['V0'])
slope1 = d1['pvals']['b']
slope_err1 = d1['perrs']['b']
offset1 = np.mean(v01[63:125])
means1 = np.mean(v01[178:184])-offset1, np.mean(v01[189:196])-offset1, np.mean(v01[203:209])-offset1, np.mean(v01[221:228])-offset1, np.mean(v01[236:243])-offset1, np.mean(v01[253:258])-offset1, np.mean(v01[267:275])-offset1
#means = np.array(np.mean(ploc[])
stdevs1 = np.std(v01[178:184]), np.std(v01[189:196]), np.std(v01[203:209]), np.std(v01[221:228]), np.std(v01[236:243]), np.std(v01[253:258]), np.std(v01[267:275])
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']
conc1 = np.array([0.0722, 0.3560, 0.7052, 3.460, 6.846, 33.57, 66.47])
# fig, ax = plt.subplots()
# x_vals = np.array([100, 500, 1000, 5000, 10000, 50000, 100000])

#trial 2 PS
fn2 = 'DOPS/20190724_003.dat'
d2 = json.load(open(fn2, "r"))
t2 = np.array([d2['data'][j]['time'] for j in range(len(d2['data']))])
pt2 = np.arange(len(d2['data']))
v02 = np.array(d2['pvals']['V0'])
v0err2 = np.array(d2['perrs']['V0'])
slope2 = d2['pvals']['b']
slope_err2 = d2['perrs']['b']
offset2 = np.mean(v02[96:122])
means2 = np.mean(v02[164:168])-offset2, np.mean(v02[175:181])-offset2, np.mean(v02[189:195])-offset2, np.mean(v02[208:215])-offset2, np.mean(v02[224:234])-offset2, np.mean(v02[240:248])-offset2, np.mean(v02[250:259])-offset2
#means = np.array(np.mean(ploc[])
stdevs2 = np.std(v02[164:168]), np.std(v02[175:181]), np.std(v02[189:195]), np.std(v02[208:215]), np.std(v02[224:234]), np.std(v02[240:248]), np.std(v02[250:259])
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']
conc2 = np.array([0.06930, 0.3428, 0.6801, 3.343, 6.617, 32.47, 64.07])
# fig, ax = plt.subplots()
# x_vals = np.array([100, 500, 1000, 5000, 10000, 50000, 100000])

#trial 3 PS
fn3 = 'DOPS/20190724_005.dat'
d3 = json.load(open(fn3, "r"))
t3 = np.array([d3['data'][j]['time'] for j in range(len(d3['data']))])
pt3 = np.arange(len(d3['data']))
v03 = np.array(d3['pvals']['V0'])
v0err3 = np.array(d3['perrs']['V0'])
slope3 = d3['pvals']['b']
slope_err3 = d3['perrs']['b']
offset3 = np.mean(v03[81:93])
means3 = np.mean(v03[87:98])-offset3, np.mean(v03[105:110])-offset3, np.mean(v03[119:127])-offset3, np.mean(v03[133:140])-offset3, np.mean(v03[146:156])-offset3, np.mean(v03[161:171])-offset3, np.mean(v03[172:181])-offset3
#means = np.array(np.mean(ploc[])
stdevs3 = np.std(v03[87:98]), np.std(v03[105:110]), np.std(v03[119:127]), np.std(v03[133:140]), np.std(v03[146:156]), np.std(v03[161:171]), np.std(v03[172:181])
# you can similarly get the slopes using d['pvals']['b'] and d['perrs']['b']
conc3 = np.array([0.07647, 0.3778, 0.7475, 3.662, 7.240, 35.46, 70.12])
# fig, ax = plt.subplots()
# x_vals = np.array([100, 500, 1000, 5000, 10000, 50000, 100000])

np.savetxt('dops1.txt', np.vstack((conc1, means1, stdevs1)).T)
np.savetxt('dops2.txt', np.vstack((conc2, means2, stdevs2)).T)
np.savetxt('dops3.txt', np.vstack((conc3, means3, stdevs3)).T)

#graphs
plt.errorbar(conc1, means1, stdevs1, fmt='d')
plt.errorbar(conc2, means2, stdevs2, fmt='d')
plt.errorbar(conc3, means3, stdevs3, fmt='d')
plt.xscale('log')
plt.xlabel('[CaCl2] (mM)')
plt.ylabel('Voltage Offset (mV)')
#labels = [item.get_text() for item in ax.get_xticklabels()]
#labels[1] = '20 mM'

#ax.set_xticklabels(labels)

plt.show()