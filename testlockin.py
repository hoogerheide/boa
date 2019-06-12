import harmInstrI
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