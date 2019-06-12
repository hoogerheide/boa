import harmInstrI_SR860_update as harmInstrI
import numpy
lockin = harmInstrI.SR860()
# print lockin.dcVoltage
lockin.dcVoltage = 0.2
print(lockin.ctrl.query("AUXV? 1"))
print(harmInstrI.str2num(lockin.ctrl.query("FREQ?"))[0])
# SR860 scan commands, SCNDC, SCNEND 2 (xii) with data capture commmands (xiii) for DC voltage offset
        # Capture rate- one-shot capture, up-down scan until scan is does (then stop) OR number of samples
        # Read in capture points, CAPTUREELEN, max capture rate 2441.41 Hz (less than possibe max 9.3 kHz)
        # Overall length >4Mb, 55s length 2441.41.*55 total points, buffer length 1074 kbytes
        # Capturestop command- stop when it becomes full, Capturelen to 1074 kbytes, capturecfg xy, 
        # captureratemax, actual rate /n^2, fill buffer highest sampling rate over 55s
        # Captureget - read data back
        # Could structure this as a function w/ dictionary of values

# Set parameter values
scnTime = 15
# Seconds- could make function to convert from minutes
scnStart = 0.05
scnEnd = -0.05
# Voltage- -5.00V < V < 5.00V
scnInt = 0
# Seconds/msec- 0 = 8ms 


lockin.ctrl.write("SCNRST")
# Resets scan regardless of state, resets to begin parameter (SCNENBL may be redundant)
# lockin.ctrl.write("SCNENBL OFF")
lockin.ctrl.write("SCNPAR REFD")
# Set scan parameter to REFDc (reference DC)
lockin.ctrl.write("SCNLOG LIN")
# Set scan type to linear
lockin.ctrl.write("SCNEND UP")
# Set scan end mode to UP (updown)
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
lockin.ctrl.write("SCNRUN")
print(lockin.ctrl.write("SCNSTATE?"))

""" # SR860 capture commands
lockin.ctrl.write("CAPTURECFG XY")
# Set capture configuration to X and Y
lockin.ctrl.write("CAPTURERATEMAX 325000")
# Set capture configuration to max rate in Hz
    # lockin.ctrl.write("CAPTURERATEMAX" + `maxRate`)
    # maxRate = 325000
lockin.ctrl.write("CAPTURERATE 0")
# Set capture rate to maximum rate /2^n for n=0
    # lockin.ctrl.write("CAPTURERATE" + `n`)
    # n = 0
lockin.ctrl.write("CAPTURESTART ONE, IMM")
# Set data capture for oneshot (fills buffer once) acquisition with immediate start (could implement a hardware trigger)
print(lockin.ctrl.query("CAPTURERATEMAX?"))
print(lockin.ctrl.query("CAPTURELEN?"))"""
