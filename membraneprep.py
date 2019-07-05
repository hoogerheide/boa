""" This is a script to prime the membrane for formation """


import harmInstrI_SR860_update as harmInstrI
lockin = harmInstrI.SR860()
lockin.sensitivity = "100 mV/nA"
lockin.harmonic = 1
lockin.close()