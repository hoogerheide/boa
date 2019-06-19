import harmInstrI_SR860_update as harmInstrI
import numpy as np
import time
import struct
import matplotlib
matplotlib.rcParams['backend'] = 'wxagg'
import matplotlib.pyplot as plt

#  This is a function to format the query output from string--> list of strings --> list of numbers
def str2num(str):
    li = str[:-1].split(",")
    return map(float,li)

lockin = harmInstrI.SR860()
# print lockin.dcVoltage
lockin.dcVoltage = 0.2
# print(lockin.ctrl.query("AUXV? 1"))
# print(harmInstrI.str2num(lockin.ctrl.query("FREQ?"))[0])

# SR860 scan commands, SCNDC, SCNEND 2 (xii) with data capture commmands (xiii) for DC voltage offset
        # Capture rate- one-shot capture, up-down scan until scan is does (then stop) OR number of samples
        # Read in capture points, CAPTUREELEN, max capture rate 2441.41 Hz (less than possibe max 9.3 kHz)
        # Overall length >4Mb, 55s length 2441.41.*55 total points, buffer length 1074 kbytes
        # Capturestop command- stop when it becomes full, Capturelen to 1074 kbytes, capturecfg xy, 
        # captureratemax, actual rate /n^2, fill buffer highest sampling rate over 55s
        # Captureget - read data back
        # Could structure this as a function w/ dictionary of values

# Calculations for scan parameters
        # IEEE 488 parameter 1.8MB/s, sampling cutoff frequency guesstimate ~50Hz
# Set parameter values
# Calculate capture length based on scan time
scnTime = 15

# Seconds- could make function to convert from minutes
scnStart = 0.05
scnEnd = -0.05
# Voltage- -5.00V < V < 5.00V
scnInt = 0
# Seconds or msec- 0 = 8ms 

filterdict = {"1 us": (0, 1e-6),
                           "3 us": (1, 3e-6),
                           "10 us": (2, 10e-6),
                           "30 us": (3, 30e-6),
                           "100 us": (4, 100e-6),
                           "300 us": (5, 300e-6),
                           "1 ms": (6, 1e-3),
                           "3 ms": (7, 3e-3),
                           "10 ms": (8, 10e-3),
                           "30 ms": (9, 30e-3),
                           "100 ms": (10, 100e-3),
                           "300 ms": (11, 300e-3),
                           "1 s": (12, 1.),
                           "3 s": (13, 3.),
                           "10 s": (14, 10.),
                           "30 s": (15, 30.),
                           "100 s": (16, 100),
                           "300 s": (17, 300),
                           "1 ks": (18, 1e3),
                           "3 ks": (19, 3e3),
                           "10 ks": (20, 10e3),
                           "30 ks": (21, 30e3)}
tConstant=lockin.filterdict[lockin.inv_filterdict[int(lockin.ctrl.query("OFLT?")[:-1])]][1]

maxRate = str2num(lockin.ctrl.query("CAPTURERATEMAX?"))[0]
# Set capture rate to max rate in Hz, where T=100ms, maxcapture=325kHz
# 1uS < T < 10uS, maxcapture=1250000
# Could create dictionary to select- write function to calculate capLength from maxRate
capdict = {"1 to 10 us": (0, 1.25e6),
                           "30 us": (1, 625e3),
                           "100 us": (2, 325e3),
                           "300 us": (3, 156.25e3),
                           "1 ms": (4, 78.125e3),
                           "3 to 10 ms": (5, 39.0625e3),
                           "30 ms": (6, 9765.62),
                           "100 ms": (7, 2441.41),
                           "300 ms": (8, 1220.7),
                           "1 s": (9, 305.18),
                           "3 s to 30 ks": (10, 152.59)}
fCuttoff= 5/tConstant
# Write function to calculate fCuttoff for low-pass filter based on tConstant
maxArray = maxRate/2**np.arange(21)
nFactor = np.where(maxArray>fCuttoff)[0][-1]
capLength = np.ceil(maxArray[nFactor]*scnTime*8/1000)


lockin.ctrl.write("SCNRST")
# Resets scan regardless of state, resets to begin parameter (SCNENBL may be redundant)
# lockin.ctrl.write("SCNENBL OFF")
lockin.ctrl.write("SCNPAR REFD")
# Set scan parameter to REFDc (reference DC)
lockin.ctrl.write("SCNLOG LIN")
# Set scan type to linear
lockin.ctrl.write("SCNEND 0")
# Set scan end mode to UP (updown), RE (repeat), ON (once)
lockin.ctrl.write("SCNSEC " + `scnTime`)
# Set scan time to x seconds (scnTime)
lockin.ctrl.write("SCNDCATTN 0")
# Set dc output attenuator mode to auto 0 or fixed 1
lockin.ctrl.write("SCNDC BEG, " + `scnStart`)
lockin.ctrl.write("SCNDC END, " + `scnEnd`)
# Set beginning (BEG) and end (END) dc reference amplitude to V, where -5.00V < V < 5.00V
lockin.ctrl.write("SCNINRVL " + `scnInt`)
# Set parameter update interval 0 <= scnInt <= 16 according to numeric table (manual pg 129)
lockin.ctrl.write("SCNENBL ON")
# Set scan parameter to begin value but does not start scan

# SR860 capture commands
lockin.ctrl.write("CAPTURECFG XY")
# Set capture configuration to X and Y
#lockin.ctrl.write("CAPTURERATEMAX " + `maxRate`)
# Set capture configuration to max rate in Hz
lockin.ctrl.write("CAPTURERATE " + `nFactor`)
# Set capture rate to maximum rate /2^n for n=0
lockin.ctrl.write("CAPTURELEN " + `capLength`)
print(lockin.ctrl.query("CAPTURERATEMAX?"))
print(lockin.ctrl.query("CAPTURELEN?"))
print(lockin.ctrl.query("CAPTUREBYTES?"))

starttime = time.time()
lockin.ctrl.write("SCNRUN")
lockin.ctrl.write("CAPTURESTART ONE, IMM")

# Begin data capture only after SCNSTATE reflects done
state = 0
while state < 4:
    state=str2num(lockin.ctrl.query("SCNSTATE?"))[0]
    print(state)

lockin.ctrl.write("CAPTURESTOP")
# Set data capture for oneshot (fills buffer once) acquisition with immediate start (could implement a hardware trigger)

# Data capture results
# Unpacking binary data- time start now

lockin.ctrl.write("CAPTUREGET? 0,%i" % capLength)
t = time.time()-starttime
buf = lockin.ctrl.read_raw()
hdr = struct.unpack_from('<cc', buf)
datalength = struct.unpack_from('<' + 'c'*int(hdr[1]), buf, 2)
#print(hdr, '<%if' % (int("".join(datalength))/4))
data = struct.unpack_from('<%if' % (int("".join(datalength))/4), buf, 2 + int(hdr[1]))
Y = np.array(data[1::2])
X = np.array(data[0::2])
# Chop off trailing data in buffer after stop capture
idx = np.where(t>scnTime)[0]
X=X[0:idx]
Y=Y[0:idx]
tStep = scnTime/idx
vTime = np.arange(0,scnTime,tStep)
vStep = np.arange(scnStart,scnEnd,tStep)

# Plot 
"""plt.plot(np.arange(len(X))/maxArray[nFactor], X)
plt.plot(np.arange(len(Y))/maxArray[nFactor], Y)
plt.show()"""