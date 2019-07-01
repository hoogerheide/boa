""" Library of instrument objects for second harmonic analysis. Only the SR830
    is used in this configuration. Ignore the filter.
    
    TODO: create an initialization routine for the SR830 and filter that will
          take them from an arbitrary state and make them ready for second
          harmonic analysis."""

import numpy
from visa import ResourceManager
import time
import struct

#  This is a function to format the query output from string--> list of strings --> list of numbers
def str2num(str):
    li = str[:-1].split(",")
    return map(float,li)

rm=ResourceManager()

class SR830(object):
    """ Stanford Research Systems SR830 instrument class 
        Manual: https://www.thinksrs.com/downloads/pdfs/manuals/SR830m.pdf """
    
    # Define dictionary values for SR830 object

    sensitivitydict =  {"2 nV/fA": 0,
                                "5 nV/fA": 1,
                                "10 nV/fA": 2,
                                "20 nV/fA": 3,
                                "50 nV/fA": 4,
                                "100 nV/fA": 5,
                                "200 nV/fA": 6,
                                "500 nV/fA": 7,
                                "1 uV/pA": 8,
                                "2 uV/pA": 9,
                                "5 uV/pA": 10,
                                "10 uV/pA": 11,
                                "20 uV/pA": 12,
                                "50 uV/pA": 13,
                                "100 uV/pA": 14,
                                "200 uV/pA": 15,
                                "500 uV/pA": 16,
                                "1 mV/nA": 17,
                                "2 mV/nA": 18,
                                "5 mV/nA": 19,
                                "10 mV/nA": 20,
                                "20 mV/nA": 21,
                                "50 mV/nA": 22,
                                "100 mV/nA": 23,
                                "200 mV/nA": 24,
                                "500 mV/nA": 25,
                                "1 V/uA": 26}
                                
    filterdict = {"10 us": (0, 10e-6),
                            "30 us": (1, 30e-6),
                            "100 us": (2, 100e-6),
                            "300 us": (3, 300e-6),
                            "1 ms": (4, 1e-3),
                            "3 ms": (5, 3e-3),
                            "10 ms": (6, 10e-3),
                            "30 ms": (7, 30e-3),
                            "100 ms": (8, 100e-3),
                            "300 ms": (9, 300e-3),
                            "1 s": (10, 1.),
                            "3 s": (11, 3.),
                            "10 s": (12, 10.),
                            "30 s": (13, 30.),
                            "100 s": (14, 100.),
                            "300 s": (15, 300.),
                            "1 ks": (16, 1e3),
                            "3 ks": (17, 3e3),
                            "10 ks": (18, 10e3),
                            "30 ks": (19, 30e3)}
                            
    filterslopedict = {"6 dB/oct": 0,
                                    "12 dB/oct": 1,
                                    "18 dB/oct": 2,
                                    "24 dB/oct": 3}
                                    
    reservedict = {"high": 0,
                                    "normal": 1,
                                    "low noise": 2}
    
    def __init__( self, address=8 ):
        self.ctrl = rm.open_resource( "GPIB::%s" % address )
        self.ctrl.open()
        
        self.ctrl.write("OUTX 1") # set output to GPIB
        
        # Set dictionaries of values
        self.gains = dict({0.5: 0.05, 1: 0.1, 1.5: 0.2, 2.0: 0.5, 2.5: 1, 3.0: 2, 3.5: 5, 4.0: 10, 4.5: 20, 5.0: 50, 5.5: 100, 6.0: 200, 6.5: 500})
        
        """self.filterdict = filterdict
        self.filterslopedict = filterslopedict
        self.sensitivitydict = sensitivitydict
        self.reservedict = reservedict"""
        
        self.inv_filterdict = {v[0]: k for k, v in self.filterdict.items()}
        self.inv_filterslopedict = {v: k for k, v in self.filterslopedict.items()}
        self.inv_sensitivitydict = {v: k for k, v in self.sensitivitydict.items()}
        self.inv_reservedict = {v: k for k, v in self.reservedict.items()}   
        
       # self.ReadValues()
        
    def Initialize(self):
        # Initialize lockin to correct state for second harmonics measurement
        self.ctrl.write("RSLP 0") # set reference to "SINE"
        self.ctrl.write("FMOD 1") # set reference to internal
        self.ctrl.write("ISRC 2") # set to perform current measurement at I*10^6 sensitivity

    def ReadValues(self):
        # Populate properties with values read from instrument
        # Update: .ask replaced with .query, str2num function to format output
        self._filter = self.inv_filterdict[int(self.ctrl.query("OFLT?")[:-1])]
        self._filterslope = self.inv_filterslopedict[int(self.ctrl.query("OFSL?")[:-1])]
        self._sensitivity = self.inv_sensitivitydict[int(self.ctrl.query("SENS?")[:-1])]
        self._reserve = self.inv_reservedict[int(self.ctrl.query("RMOD?")[:-1])]
        self._frequency = str2num(self.ctrl.query("FREQ?")[0])
            # alternatively- str2num(self._frequency)
        self._amplitude = str2num(self.ctrl.query("SLVL?")[0])
            # alternatively- str2num(self._amplitude)
        self._harmonic = str2num(self.ctrl.query("HARM?")[0])
            # alternatively- str2num(self._harmonic)
        self._dcVoltage = str2num(self.ctrl.query("AUXV? 1")[0])
            # alternatively- str2num(self.dcVoltage)

    def close(self):
        """ closes the VISA instance (I think) """
        self.ctrl.close()
    
    @property
    def filter(self):
        return self._filter
        
    @filter.setter
    def filter(self, key):
        if key in self.filterdict.keys():
            self._filter = key
            self.ctrl.write("OFLT %i" % self.filterdict[key][0])
        else:
            print "Warning: filter setting %s not valid" % key
            
    @property
    def filterslope(self):
        return self._filterslope
        
    @filterslope.setter
    def filterslope(self, key):
        if key in self.filterslopedict.keys():
            self._filterslope = key
            self.ctrl.write("OFSL %i" % self.filterslopedict[key])
        else:
            print "Warning: filter slope setting %s not valid" % key
            
    @property
    def reserve(self):
        return self._reserve
        
    @reserve.setter
    def reserve(self, key):
        if key in self.reservedict.keys():
            self._reserve = key
            self.ctrl.write("RMOD %i" % self.reservedict[key])
        elif key == 'auto':
            self.autoReserve()
        else:
            print "Warning: reserve setting %s not valid" % key            
            
    @property
    def sensitivity(self):
        return self._sensitivity
        
    @sensitivity.setter
    def sensitivity(self, key):
        if key in self.sensitivitydict.keys():
            self._sensitivity = key
            self.ctrl.write("SENS %i" % self.sensitivitydict[key])
        else:
            print "Warning: sensitivity setting %s not valid" % key
                
    @property
    def frequency(self):
        return self._frequency
    
    @frequency.setter
    def frequency(self, value):
        self._frequency = value        
        self.ctrl.write("FREQ " + `value`)
    
    @property
    def amplitude(self):
        return self._amplitude
    
    @amplitude.setter
    def amplitude(self, value):
        if value > 5.:
            print "Warning: Lock-in amplitude cannot be greater than 5 V. Reducing to 5 V."
            value = 5.
            
        self._amplitude = value
        self.ctrl.write("SLVL " + `value`)
    
    @property
    def harmonic(self):
        return self._harmonic
    
    @harmonic.setter
    def harmonic(self, value):
        self._harmonic = value
        self.ctrl.write("HARM " + `value`)
    
    @property
    def dcVoltage(self):
        return self._dcVoltage
    
    @dcVoltage.setter
    def dcVoltage(self, value):
        # TODO: can't be greater than +/-10.5 V. Check for this.
        self._dcVoltage = value
        self.ctrl.write("AUXV 1, " + `value`)
        
    @property
    def phase(self):
        return self._phase
        
    @phase.setter
    def phase(self, value):
        # TODO: can't be greater than +730/-360. Check for this.
        self._phase = value
        self.ctrl.write("PHAS " + `value`)
        
    def getPhase(self):
        return str2num(self.ctrl.query("PHAS?")[0])
        
    def autoPhase(self):
        self.ctrl.write("APHS")
        
    def autoReserve(self):
        self.ctrl.write("ARSV")
        
    def hasOverload(self):
        """ Looks for an input or amplifier overload. """
        return  str2num(bool(self.ctrl.query("LIAE? 0")[0]))
    
    def gain(self, auxchannel=1):
        # not used with current measurement
        pass
        return str2num(float(self.gains[round(self.ctrl.query("OAUX? " + `auxchannel`)[0]*2)/2]))

    def measureRTheta(self):
        return str2num(self.ctrl.query('SNAP? 3,4'))
        
    def measureXY(self):
        return str2num(self.ctrl.query('SNAP? 1,2'))
        
    def sweep(self, propertyName, vals, equilInterval, numPoints, pointInterval):
        """ Attempt to sweep a given property (e.g. dcVoltage or acAmplitude).
            For second harmonic analysis, this would allow the use of the same 
            function for both second and third harmonic determination. Not
            currently used. """
                
        ival = self.__getattribute__(propertyName)
        
        data = dict()
        data['vals'] = vals
        data['X'] = numpy.empty(vals.size)
        data['Xerr'] = numpy.empty(vals.size)
        data['Y'] = numpy.empty(vals.size)
        data['Yerr'] = numpy.empty(vals.size)
        
        for i,val in enumerate(vals):
            self.__setattr__(propertyName, val)
            time.sleep(equilInterval)
            data['X'][i], data['Y'][i], data['Xerr'][i], data['Yerr'][i] = self.collectPoints(numPoints, pointInterval)
        
        self.__setattr__(propertyName, ival)
        
        return data
    
    def collectPointsXY(self, numPoints, pointInterval):
        Xs = numpy.empty(numPoints)
        Ys = numpy.empty(numPoints)
        for j in range(numPoints):
            Xs[j], Ys[j] = self.measureXY()
            time.sleep(pointInterval)
        
        return numpy.mean(Xs), numpy.mean(Ys), numpy.std(Xs), numpy.std(Ys)
        
    def collectPointsRTheta(self, numPoints, pointInterval):
        Rs = numpy.empty(numPoints)
        Ts = numpy.empty(numPoints)
        for j in range(numPoints):
            Rs[j], Ts[j] = self.measureRTheta()
            time.sleep(pointInterval)

        Ts = (Ts + 360) % 360.
        
        return numpy.mean(Rs), numpy.mean(Ts), numpy.std(Rs), numpy.std(Ts)

    def collectPoints(self, numPoints, pointInterval):
        """ For back compatibility """
        return self.collectPointsRTheta(numPoints, pointInterval)


class SR860(object):
    """ Stanford Research Systems SR860 instrument class 
        Manual: https://www.thinksrs.com/downloads/pdfs/manuals/SR860m.pdf """
    # Define dictionary values for SR860 object

    sensitivitydict =  {"1 V[uA]": 0,
                                "500 mV/nA": 1,
                                "200 mV/nA": 2,
                                "100 mV/nA": 3,
                                "50 mV/nA": 4,
                                "20 mV/nA": 5,
                                "10 mV/nA": 6,
                                "5 mV/nA": 7,
                                "2 mV/nA": 8,
                                "1 mV/nA": 9,
                                "500 uV/pA": 10,
                                "200 uV/pA": 11,
                                "100 uV/pA": 12,
                                "50 uV/pA": 13,
                                "20 uV/pA": 14,
                                "10 uV/pA": 15,
                                "5 uV/pA": 16,
                                "2 uV/pA": 17,
                                "1 uV/pA": 18,
                                "500 nV/fA": 19,
                                "200 mV/fA": 20,
                                "100 nV/fA": 21,
                                "50 nV/fA": 22,
                                "20 nV/fA": 23,
                                "10 nV/fA": 24,
                                "5 nV/fA": 25,
                                "2 nV/fA": 26,
                                "1 nV/fA": 27}
                                
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
                            
    filterslopedict = {"6 dB/oct": 0,
                                    "12 dB/oct": 1,
                                    "18 dB/oct": 2,
                                    "24 dB/oct": 3}
                                    
    reservedict = {"high": 0,
                                    "normal": 1,
                                    "low noise": 2}

    def __init__( self, address=4 ):
        self.ctrl = rm.open_resource( "GPIB::%s" % address )
        self.ctrl.open()
        # Check OUTX 1 channel address for SR860
        # self.ctrl.write("OUTX 1") # set output to GPIB - Obsolete for SR860
        # Set dictionaries of values- may need to update values
        self.gains = dict({0.5: 0.05, 1: 0.1, 1.5: 0.2, 2.0: 0.5, 2.5: 1, 3.0: 2, 3.5: 5, 4.0: 10, 4.5: 20, 5.0: 50, 5.5: 100, 6.0: 200, 6.5: 500})
        
        """self.filterdict = filterdict
        self.filterslopedict = filterslopedict
        self.sensitivitydict = sensitivitydict
        self.reservedict = reservedict"""
        
        self.inv_filterdict = {v[0]: k for k, v in self.filterdict.items()}
        self.inv_filterslopedict = {v: k for k, v in self.filterslopedict.items()}
        self.inv_sensitivitydict = {v: k for k, v in self.sensitivitydict.items()}
        self.inv_reservedict = {v: k for k, v in self.reservedict.items()}
        self.ReadValues()
        
    def Initialize(self):
        # Initialize lockin to correct state for second harmonics measurement
        self.ctrl.write("RTRG SIN") # set reference to "SINE"- double check correct command, it may be default setting
        self.ctrl.write("RSRC INT") # set reference to internal, may be obsolete with FREQINT
        self.ctrl.write("ICUR 1MEG") # set to perform current measurement at I*10^6 (Mohm) sensitivity- confirm correct state for SR860
       
    def ReadValues(self):
        # Populate properties with values read from instrument
        # Update: .ask replaced with .query, str2num function to format output
        self._filter = self.inv_filterdict[int(self.ctrl.query("OFLT?")[:-1])]
        self._filterslope = self.inv_filterslopedict[int(self.ctrl.query("OFSL?")[:-1])]
        self._sensitivity = self.inv_sensitivitydict[int(self.ctrl.query("SCAL?")[:-1])]
        # self._reserve = self.inv_reservedict[int(self.ctrl.query("RMOD?")[:-1])] - no dynamic reserve setting on SR860
        self._frequency = str2num(self.ctrl.query("FREQINT?"))[0]
        self._amplitude = str2num(self.ctrl.query("SLVL?"))[0]
        self._harmonic = str2num(self.ctrl.query("HARM?"))[0]
        self._dcVoltage = str2num(self.ctrl.query("SOFF?"))[0]
     
    def close(self):
        """ closes the VISA instance (I think) """
        self.ctrl.close()
    
    @property
    def filter(self):
        return self._filter
        
    @filter.setter
    def filter(self, key):
        if key in self.filterdict.keys():
            self._filter = key
            self.ctrl.write("OFLT %i" % self.filterdict[key][0])
        else:
            print "Warning: filter setting %s not valid" % key
            
    @property
    def filterslope(self):
        return self._filterslope
        
    @filterslope.setter
    def filterslope(self, key):
        if key in self.filterslopedict.keys():
            self._filterslope = key
            self.ctrl.write("OFSL %i" % self.filterslopedict[key])
        else:
            print "Warning: filter slope setting %s not valid" % key
            
    """ @property
    def reserve(self):
        return self._reserve """
        
    """ @reserve.setter # may be obsolete SR860
    def reserve(self, key):
        if key in self.reservedict.keys():
            self._reserve = key
            # self.ctrl.write("RMOD %i" % self.reservedict[key])
            # No dynamic reserve setting on SR860- see pg 48-49
        elif key == 'auto':
            self.autoReserve() # change to autoScale (line 462)- no auto reserve key on SR860
        else:
            print "Warning: reserve setting %s not valid" % key """            
            
    @property
    def sensitivity(self):
        return self._sensitivity
        
    @sensitivity.setter
    def sensitivity(self, key):
        if key in self.sensitivitydict.keys():
            self._sensitivity = key
            self.ctrl.write("SCAL %i" % self.sensitivitydict[key])
        else:
            print "Warning: sensitivity setting %s not valid" % key
                
    @property
    def frequency(self):
        return self._frequency
    
    @frequency.setter
    def frequency(self, value):
        self._frequency = value        
        self.ctrl.write("FREQ " + `value`)
    
    @property
    def amplitude(self):
        return self._amplitude
    
    @amplitude.setter
    def amplitude(self, value):
        if value > 5.:
            print "Warning: Lock-in amplitude cannot be greater than 5 V. Reducing to 5 V."
            value = 5.
            
        self._amplitude = value
        self.ctrl.write("SLVL " + `value`)
    
    @property
    def harmonic(self):
        return self._harmonic
    
    @harmonic.setter
    def harmonic(self, value):
        self._harmonic = value
        self.ctrl.write("HARM " + `value`)
    
    @property
    def dcVoltage(self):
        return self._dcVoltage
        # SR860 scan commands (xii) with data capture commmands (xiii) for DC voltage offset
    
    @dcVoltage.setter
    def dcVoltage(self, value):
        # TODO: can't be greater than +/-10.5 V. Check for this.
        self._dcVoltage = value
        self.ctrl.write("SOFF " + `value`) # DC output

    @property
    def phase(self):
        return self._phase
        
    @phase.setter
    def phase(self, value):
        # TODO: can't be greater than +730/-360. Check for this.
        self._phase = value
        self.ctrl.write("PHAS " + `value`)
        
    def getPhase(self):
        return str2num(self.ctrl.query("PHAS?"))[0]
        
    def autoPhase(self):
        self.ctrl.write("APHS")
        
    def autoScale(self):
        self.ctrl.write("ASCL")
        
    def hasOverload(self):
        """ Looks for an input or amplifier overload. """
        return  str2num(bool(self.ctrl.query("LIAE? 0")))[0]
    
    def gain(self, auxchannel=1):
        # not used with current measurement
        pass
        return str2num(float(self.gains[round(self.ctrl.query("OAUX? " + `auxchannel`)[0]*2)/2]))
        # OAUX 1- change to use DC offset on input, SR860 scan commands (xii) with data capture commmands (xiii) for DC voltage offset

    # measureXY , collectpointsXY functions for compatibility with GUI
    def measureXY(self):
        return str2num(self.ctrl.query('SNAP? 0,1'))
    
    def measureRTheta(self):
        return str2num(self.ctrl.query('SNAP? 2,3'))

    def collectPointsXY(self, numPoints, pointInterval):
        Xs = numpy.empty(numPoints)
        Ys = numpy.empty(numPoints)
        for j in range(numPoints):
            Xs[j], Ys[j] = self.measureXY()
            time.sleep(pointInterval)
        
        return numpy.mean(Xs), numpy.mean(Ys), numpy.std(Xs), numpy.std(Ys)

    # Function to run scan and capture, unpack and return X, Y, V
    def measureXYV(self,scnStart,scnEnd,scnTime):
        # Initialize lockin to correct state for scan measurement
        # Calculate capture length based on scan time, move to params
        # Voltage- -5.00V < V < 5.00V
        #print scnTime, type(scnTime)
        scnInt = 0
        # Seconds or msec- 0 = 8ms 
        tConstant=self.filterdict[self.inv_filterdict[int(self.ctrl.query("OFLT?")[:-1])]][1]
        #print 'tConstant', tConstant
        maxRate = str2num(self.ctrl.query("CAPTURERATEMAX?"))[0]
        #print 'maxRate', maxRate
        fCuttoff= 5/tConstant
        # Write function to calculate fCuttoff for low-pass filter based on tConstant
        maxArray = maxRate/2**numpy.arange(21)
        self.ctrl.write("SCNRST") # sets scan regardless of state, resets to begin parameter (SCNENBL may be redundant)
        self.ctrl.write("SCNPAR REFD") # Set scan parameter to REFDc (reference DC)
        self.ctrl.write("SCNLOG LIN") # Set scan type to linear
        self.ctrl.write("SCNEND 0") # Set scan end mode to UP (updown), RE (repeat), ON (once)
        self.ctrl.write("SCNSEC " + `scnTime`) # Set scan time to x seconds (scnTime)
        self.ctrl.write("SCNDCATTN 0") # Set dc output attenuator mode to auto 0 or fixed 1
        self.ctrl.write("SCNDC BEG, " + `scnStart`) 
        self.ctrl.write("SCNDC END, " + `scnEnd`) # Set beginning (BEG) and end (END) dc reference amplitude to V, where -5.00V < V < 5.00V
        self.ctrl.write("SCNINRVL " + `scnInt`) # Set parameter update interval 0 <= scnInt <= 16 according to numeric table (manual pg 129)
        self.ctrl.write("SCNENBL ON") # Set scan parameter to begin value but does not start scan

        # Initialize lockin to correct state for capture measurement
        nFactor = numpy.where(maxArray>fCuttoff)[0][-1]
        #print 'nFactor', nFactor
        capLength = numpy.ceil(maxArray[nFactor]*scnTime*8/1000)
        #print 'capLength', capLength
        self.ctrl.write("CAPTURECFG XY") # Set capture configuration to X and Y
        self.ctrl.write("CAPTURERATE " + `nFactor`) # Set capture rate to maximum rate /2^n for n=0
        self.ctrl.write("CAPTURELEN " + `capLength`) # Set capture length according to formula in params
        starttime = time.time()
        self.ctrl.write("SCNRUN")
        self.ctrl.write("CAPTURESTART ONE, IMM")

        # Begin data capture get only after SCNSTATE reflects done
        state = 0
        while state < 4:
            state=str2num(self.ctrl.query("SCNSTATE?"))[0]

        # Stop data capture before capture get commands
        self.ctrl.write("CAPTURESTOP")

        # Data capture results
        # Unpacking binary data- time start now
        t = time.time()-starttime
        self.ctrl.write("CAPTUREGET? 0,%i" % capLength)
        buf = self.ctrl.read_raw() # Read buffer contents
        #print len(buf), maxArray[nFactor], scnTime, maxArray[nFactor]*scnTime
        hdr = struct.unpack_from('<cc', buf) # Read binary header in little endian format
        datalength = struct.unpack_from('<' + 'c'*int(hdr[1]), buf, 2)
        data = struct.unpack_from('<%if' % (int("".join(datalength))/4), buf, 2 + int(hdr[1]))
        Y = numpy.array(data[1::2])
        X = numpy.array(data[0::2])
        # Chop off trailing data in buffer after stop capture
        #idx = numpy.where(t>scnTime)[0]
        idx = int(numpy.floor(maxArray[nFactor]*scnTime))
        X=X[0:idx]
        Y=Y[0:idx]
        #tStep = scnTime/idx
        #vTime = numpy.arange(0,scnTime,tStep)
        V = numpy.arange(scnStart,scnEnd,(scnEnd-scnStart)/float(idx))
        #print "End of collection", X, Y, V
        return X, Y, V
        
    def sweep(self, propertyName, vals, equilInterval, numPoints, pointInterval):
        """ Attempt to sweep a given property (e.g. dcVoltage or acAmplitude).
            For second harmonic analysis, this would allow the use of the same 
            function for both second and third harmonic determination. Not
            currently used. """
                
        ival = self.__getattribute__(propertyName)
        
        data = dict()
        data['vals'] = vals
        data['X'] = numpy.empty(vals.size)
        data['Xerr'] = numpy.zeros(vals.size)
        # Xerr- array of zeroes for compatibility
        data['Y'] = numpy.empty(vals.size)
        data['Yerr'] = numpy.zeros(vals.size)
        # Yerr- array of zeroes for compatibility
        
        for i,val in enumerate(vals):
            self.__setattr__(propertyName, val)
            time.sleep(equilInterval)
            data['X'][i], data['Y'][i], data['Xerr'][i], data['Yerr'][i] = self.collectPoints(numPoints, pointInterval)
        
        self.__setattr__(propertyName, ival)
        
        return data

    """def collectPointsXYV(self, numPoints, pointInterval):
        Xs = numpy.empty(numPoints)
        Ys = numpy.empty(numPoints)
        Vs = numpy.empty(numPoints)
        for j in range(numPoints):
            Xs[j], Ys[j], Vs[j] = self.measureXYV()
            time.sleep(pointInterval)
        
        return numpy.mean(Xs), numpy.mean(Ys), numpy.mean(Vs), numpy.std(Xs), numpy.std(Ys), numpy.std(Vs)"""
        
    def collectPointsRTheta(self, numPoints, pointInterval):
        Rs = numpy.empty(numPoints)
        Ts = numpy.empty(numPoints)
        for j in range(numPoints):
            Rs[j], Ts[j] = self.measureRTheta()
            time.sleep(pointInterval)

        Ts = (Ts + 360) % 360.
        
        return numpy.mean(Rs), numpy.mean(Ts), numpy.std(Rs), numpy.std(Ts)

    def collectPoints(self, numPoints, pointInterval):
        """ For back compatibility """
        return self.collectPointsRTheta(numPoints, pointInterval)

class lockinZapper(SR830):
    """ Stripped-down lockin for zapping membranes. Now obsolete. """
    sensitivitydict =  {"2 nV/fA": 0,
                                "5 nV/fA": 1,
                                "10 nV/fA": 2,
                                "20 nV/fA": 3,
                                "50 nV/fA": 4,
                                "100 nV/fA": 5,
                                "200 nV/fA": 6,
                                "500 nV/fA": 7,
                                "1 uV/pA": 8,
                                "2 uV/pA": 9,
                                "5 uV/pA": 10,
                                "10 uV/pA": 11,
                                "20 uV/pA": 12,
                                "50 uV/pA": 13,
                                "100 uV/pA": 14,
                                "200 uV/pA": 15,
                                "500 uV/pA": 16,
                                "1 mV/nA": 17,
                                "2 mV/nA": 18,
                                "5 mV/nA": 19,
                                "10 mV/nA": 20,
                                "20 mV/nA": 21,
                                "50 mV/nA": 22,
                                "100 mV/nA": 23,
                                "200 mV/nA": 24,
                                "500 mV/nA": 25,
                                "1 V/uA": 26}
                                
    filterdict = {"10 us": (0, 10e-6),
                            "30 us": (1, 30e-6),
                            "100 us": (2, 100e-6),
                            "300 us": (3, 300e-6),
                            "1 ms": (4, 1e-3),
                            "3 ms": (5, 3e-3),
                            "10 ms": (6, 10e-3),
                            "30 ms": (7, 30e-3),
                            "100 ms": (8, 100e-3),
                            "300 ms": (9, 300e-3),
                            "1 s": (10, 1.),
                            "3 s": (11, 3.),
                            "10 s": (12, 10.),
                            "30 s": (13, 30.),
                            "100 s": (14, 100.),
                            "300 s": (15, 300.),
                            "1 ks": (16, 1e3),
                            "3 ks": (17, 3e3),
                            "10 ks": (18, 10e3),
                            "30 ks": (19, 30e3)}
                            
    filterslopedict = {"6 dB/oct": 0,
                                    "12 dB/oct": 1,
                                    "18 dB/oct": 2,
                                    "24 dB/oct": 3}
                                    
    reservedict = {"high": 0,
                                    "normal": 1,
                                    "low noise": 2}
    
    def __init__( self, address=8 ):
        self.ctrl = rm.open_resource( "GPIB::%s" % address )
        self.ctrl.open()
        self.ctrl.write("OUTX 1") # set output to GPIB

        # Set dictionaries of values
        self.gains = dict({0.5: 0.05, 1: 0.1, 1.5: 0.2, 2.0: 0.5, 2.5: 1, 3.0: 2, 3.5: 5, 4.0: 10, 4.5: 20, 5.0: 50, 5.5: 100, 6.0: 200, 6.5: 500})
        
        """self.filterdict = filterdict
        self.filterslopedict = filterslopedict
        self.sensitivitydict = sensitivitydict"""
        
        self.inv_filterdict = {v[0]: k for k, v in self.filterdict.items()}
        self.inv_filterslopedict = {v: k for k, v in self.filterslopedict.items()}
        self.inv_sensitivitydict = {v: k for k, v in self.sensitivitydict.items()}
        
        # Populate properties with values read from instrument
        self._filter = self.inv_filterdict[int(self.ctrl.query("OFLT?")[:-1])]
        self._filterslope = self.inv_filterslopedict[int(self.ctrl.query("OFSL?")[:-1])]
        self._sensitivity = self.inv_sensitivitydict[int(self.ctrl.query("SENS?")[:-1])]
        self._frequency = str2num(self.ctrl.query("FREQ?")[0])
        self._amplitude = str2num(self.ctrl.query("SLVL?")[0])
        self._harmonic = str2num(self.ctrl.query("HARM?")[0])
        self._dcVoltage = str2num(self.ctrl.query("AUXV? 0")[0])
        
class FD9002(object):
    """ Frequency Devices 9002 filter instrument class """
    
    # Write-only at the moment. Possibly not robust if filter settings change
    # significantly.
    # NOT USED for current-only measurement
    # TODO: actually learn how to communicate with the filter (read values, etc.)
    
    def __init__(self, address=30):
        self.ctrl = rm.open_resource( "GPIB::%s::INSTR" % address, write_termination="\x13")
        self.ctrl.open()
        #self._frequency = None
    
    @property
    def frequency(self):
        return self._frequency
    
    @frequency.setter
    def frequency(self, value):
        self._frequency = value
        try:        
            self.ctrl.write("\x11\x0F" + `0` + "\x3B")
            time.sleep(0.05)
            self.ctrl.write("\x11\x0F" + `value` + "\x3B")
        except:
            print 'Warning: filter not set'
        
    def close(self):
        """ closes the VISA instance (I think) """
        self.ctrl.close()