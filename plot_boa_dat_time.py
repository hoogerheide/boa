import numpy as np
import json
import sys
#sys.path.append('G:\\My Drive\\Software\\boa')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

if len(sys.argv) > 1:
    
    fn = sys.argv[1]

else:
    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()

    fn = filedialog.askopenfilename()

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 1.5

def analyzescandata(x, y):
    p0X, pcovX = np.polyfit(x, y,1, full=False, cov=True)

    return np.polyval(p0X, x), p0X

d = json.load(open(fn, "r"))

Vac = d['params']['acAmplitude']*np.sqrt(2)/1000.
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
#slopex = np.array([analyzescandata(dd['V'], dd['X'])[1][0] for dd in d['data']])
#slopey = np.array([analyzescandata(dd['V'], dd['Y'])[1][0] for dd in d['data']])

fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(16,10), gridspec_kw={'height_ratios': [1,1,0.2]})
gs = axs[0,0].get_gridspec()
# remove the underlying axes
for ax in axs[:2, :2].flatten():
    ax.remove()
axbig = fig.add_subplot(gs[:2, :2])

C_ax = axs[0,2]
slope_ax = axs[0,3]
harm2_ax = axs[1,2]
harm3_ax = axs[1,3]

for ax in axs[-1,:].flatten():
    ax.remove()
ax_slider = fig.add_subplot(gs[-1, 1:3])
ax_button_left = fig.add_subplot(gs[-1, 0])
ax_button_right = fig.add_subplot(gs[-1, -1])

[linex] = axbig.plot(d['data'][0]['V'], d['data'][0]['X'], label='X', alpha=0.7, color='C1')
[linexth] = axbig.plot(d['data'][0]['V'], analyzescandata(d['data'][0]['V'], d['data'][0]['X'])[0], color=linex.get_color(), linewidth=3)
[liney] = axbig.plot(d['data'][0]['V'], d['data'][0]['Y'], label='Y', alpha=0.7, color='C2')
[lineyth] = axbig.plot(d['data'][0]['V'], analyzescandata(d['data'][0]['V'], d['data'][0]['Y'])[0], color=liney.get_color(), linewidth=3)
axbig.axvline(0, linestyle='--', color='0.1')
axbig.axhline(0, linestyle='--', color='0.1')
axbig.set_title('Time (s): %0.1f' % float(d['data'][0]['time']))
axbig.set_xlabel('Voltage (mV)')
axbig.set_ylabel('Current (pA)')
axbig.legend(loc=0)

C_ax.plot(t, capacitance)
C_ax.set_xlabel('Time (s)')
C_ax.set_ylabel('Capacitance (pF)')
[pointC] = C_ax.plot(t[0], capacitance[0], 'o', fillstyle='none', mec='red', mew=2, markersize=10)

slope_ax.errorbar(t, slope, slope_err)
#slope_ax.plot(t, slopex)
#slope_ax.plot(t, slopey)
slope_ax.set_xlabel('Time (s)')
slope_ax.set_ylabel('Slope (pA/mV)')
[pointslope] = slope_ax.plot(t[0], slope[0], 'o', fillstyle='none', mec='red', mew=2, markersize=10)

harm2_ax.errorbar(t, v0, v0err)
harm2_ax.set_xlabel('Time (s)')
harm2_ax.set_ylabel('Second harmonic offset (mV)')
[pointv0] = harm2_ax.plot(t[0], v0[0], 'o', fillstyle='none', mec='red', mew=2, markersize=10)

harm3_ax.plot(t, h3)
harm3_ax.set_xlabel('Time (s)')
harm3_ax.set_ylabel('Third harmonic amplitude (pA)')
[pointh3] = harm3_ax.plot(t[0], h3[0], 'o', fillstyle='none', mec='red', mew=2, markersize=10)

# Slider
amp_slider = Slider(ax_slider, 'Point', 0, len(d['data']), valinit=0)
# use valstep=np.arange(len(d['data'])) for Matplotlib > 3.4

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    i = int(round(val))
    linex.set_data(d['data'][i]['V'], d['data'][i]['X'])
    #linex.set_ydata(d['data'][i]['X'])
    liney.set_xdata(d['data'][i]['V'])
    liney.set_ydata(d['data'][i]['Y'])
    linexth.set_xdata(d['data'][i]['V'])
    linexth.set_ydata(analyzescandata(d['data'][i]['V'], d['data'][i]['X'])[0])
    lineyth.set_xdata(d['data'][i]['V'])
    lineyth.set_ydata(analyzescandata(d['data'][i]['V'], d['data'][i]['Y'])[0])
    pointC.set_xdata(t[i])
    pointC.set_ydata(capacitance[i])
    pointslope.set_xdata(t[i])
    pointslope.set_ydata(slope[i])
    pointv0.set_xdata(t[i])
    pointv0.set_ydata(v0[i])
    pointh3.set_xdata(t[i])
    pointh3.set_ydata(h3[i])
    axbig.set_title('Time (s): %0.1f' % float(d['data'][i]['time']))
    axbig.relim()
    axbig.autoscale(enable=True, axis='both')
    fig.canvas.draw()    

amp_slider.on_changed(sliders_on_changed)

# Add a button for changing the parameters
left_button = Button(ax_button_left, '<', color='0.975', hovercolor='gray')
right_button = Button(ax_button_right, '>', color='0.975', hovercolor='gray')
def left_button_on_clicked(mouse_event):
    amp_slider.set_val(int(round(max(amp_slider.val - 1, 0))))
    #amp_slider.set_val(amp_slider.val - 1)
def right_button_on_clicked(mouse_event):
    amp_slider.set_val(int(round(min(amp_slider.val + 1, len(d['data'])))))
left_button.on_clicked(left_button_on_clicked)
right_button.on_clicked(right_button_on_clicked)

fig.tight_layout()
plt.show()

