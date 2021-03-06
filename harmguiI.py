# -*- coding: utf-8 -*-
""" Second harmonic analysis software version 1b.0.0
    David P. Hoogerheide, 27 February 2015
    
    Updated 2 March 2015: Added "Edit tag" button
    1.0.1 (3/10/2015): Changed fitfunc from quadratic to cubic polynomial
    1.0.2 (3/10/2015): Added event handler for worker thread completion and
                        status bar updates so users don't get impatient when
                        nothing happens after hitting "stop"
    1.1.0 (3/10/2015): Incorporate harmlibinstr into worker thread class.
    1.1.1 (3/10/2015): Add lockin filter and sensitivity controls
    1.1.2 (3/11/2015): Turn txtStdOut into list control for data
    1.1.3 (5/5/2015): Changed model, added error collection for C and harm3.
    1b.0.0 (5/11/2015): Changed to be lockin-only measurement (no axopatch or filter)
    1b.1.0 (5/11/2015): Added zap button
    1b.1.1 (5/15/2015): Changed way fit is done to ensure there are no zeros in Rerr;
                        also changed default current setting to I*10^6.
    1b.1.2 (5/18/2015): Now record phase of 3rd harmonic measurement
    1b.2.0 (5/18/2015): Version that records X and Y instead of R and Theta
    1b.2.1 (5/27/2015): Fits X and Y as lines separately and calculates second harmonic
                        slope and offset from these lines
    1b.2.2 (5/28/2015): Break out analysis function to analysisXYData for use by libraries;
                        also add "Set Parameters on Lockin" button
    1b.2.3 (6/8/2015): Add reserve option. Also change harmInstrI library so init 
                        function doesn't change anything. Use Initialize() instead.
                        Also new harmInstrI.ReadValues() that reads relevant parameters
                        from the instrument.
    
    Requires files: harmInstrI.py (instrument library)
                    defaultI.json (default parameters)
"""

import wx
import wxmplot
import numpy as np
import sys
import time
import glob
from scipy.optimize import curve_fit
import json
from threading import Thread

import harmInstrI as harmInstr

class NumpyAwareJSONEncoder(json.JSONEncoder):
    """ Utility to write numpy array objects to JSON file. """
    def default(self, obj):
        if isinstance(obj, np.ndarray) and obj.ndim == 1:
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)        

def fitfunc(V, b, V0):
    """Fit function for second harmonic offset.
    Assumed to be a symmetric function around V0 with parabolic (a),
    linear (b), and constant (c) components. Technically "a" should
    be very small, but it can help deal with the complex functional
    form around V0.
    
    Note 3/10/2015: theoretically the second harmonic should only contain linear
    and cubic components. Origin of offset is also unclear.
    
    Note 5/18/2015: When using X and Y instead of R and Theta, this is a simple line.
    """
    #return a*abs(V-V0)**3 + b*abs(V-V0) + c
    
    # new function 5/5/2015
    #return np.max(np.vstack((b*abs(V-V0) + a, c + 0*abs(V-V0))), axis=0)
    
    # new function 5/18/2015
    return b*(V-V0)

def analyzeXYData(V, X, Xerr, Y, Yerr):
    """ Function to produce fit values for X, Y data. Broken out of the mainFrame
    class so it can be used in a library format.
    
    Inputs: V (array of voltages), X (array of X values), Xerr (array of X errors), Y, Yerr
    Returns: pval (slope, V0) and perr (errors in slope, V0). """
    
    # Make sure there are no zeros in error   
    sigs = np.array(Xerr)
    sigs[sigs==0] = np.mean(sigs[sigs>0])
    
    # get starting values without using error bars    
    p0X, pcovX = np.polyfit(V, X,1, full=False, cov=True)
    
    # attempt to fitting using errors. If doesn't work, revert to non-error fit
    try:
        pvalX, pcovX = curve_fit(fitfunc, V, X,
            p0=(p0X[0], -p0X[1]/p0X[0]),
            sigma=sigs,
            absolute_sigma=True,
            maxfev=2000)
    except:
        pvalX = p0X

    # Calculate estimate of errors from covariance matrix
    perrX = np.sqrt(np.diag(pcovX))
    
    # Fit Y curve
    sigs = np.array(Yerr)
    sigs[sigs==0] = np.mean(sigs[sigs>0])
    
    p0Y, pcovY = np.polyfit(V, Y, 1, full=False, cov=True)

    try:
        pvalY, pcovY = curve_fit(fitfunc, V, Y,
            p0=(p0Y[0], -p0Y[1]/p0Y[0]),
            sigma=sigs,
            absolute_sigma=True,
            maxfev=2000)
    except:
        pvalY = p0Y

    perrY = np.sqrt(np.diag(pcovY))

    # Calculate slope and offset from X, Y slopes    
    pval = np.array([np.sqrt(pvalX[0]**2 + pvalY[0]**2), (pvalX[0]**2*pvalX[1]+pvalY[0]**2*pvalY[1])/(pvalX[0]**2 + pvalY[0]**2)])
    perr = np.zeros(pval.shape)
    perr[0] = np.sqrt((2*pvalX[0]/pval[0]*perrX[0])**2+(2*pvalY[0]/pval[0]*perrY[0])**2)
    perr[1] = np.sqrt((2*pvalX[0]*(pvalX[1]-pval[1])*perrX[0]/(pvalX[0]**2 + pvalY[0]**2))**2+
                        (2*pvalY[0]*(pvalY[1]-pval[1])*perrY[0]/(pvalX[0]**2 + pvalY[0]**2))**2 +
                        (pvalX[0]**2*perrX[1]/(pvalX[0]**2 + pvalY[0]**2))**2 +
                        (pvalY[0]**2*perrY[1]/(pvalX[0]**2 + pvalY[0]**2))**2)
                        
    return pval, perr, pvalX, pvalY
    
"""
    WORKER THREAD SECTION
    The worker thread is used to collect data in the background so the gui
    stays responsive. The EVT_RESULT event and ResultEvent class are used
    to pass data to the gui whenever new data are obtained.
    """

class RedirectText(object):
    """ Define a simple class to redirect data output from worker thread """
    def __init__(self,aWxTextCtrl):
        self.out=aWxTextCtrl
 
    def write(self, string):
        wx.CallAfter(self.out.AppendText, string)

EVT_RESULT_ID = wx.NewId()
EVT_THREAD_DONE = wx.NewId()
EVT_STATUS_ID = wx.NewId()

def EVT_RESULT(win, func):
    """Define Result Event for use by parent thread"""
    win.Connect(-1, -1, EVT_RESULT_ID, func)
    
def EVT_STATUS(win, func):
    """Define Result Event for use by parent thread"""
    win.Connect(-1, -1, EVT_STATUS_ID, func)
    
def EVT_DONE(win, func):
    """Define worker completion event """
    win.Connect(-1, -1, EVT_THREAD_DONE, func)

class ResultEvent(wx.PyEvent):
    """Simple event to carry arbitrary result data."""
    def __init__(self, data):
        """Init Result Event."""
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_RESULT_ID)
        self.data = data

class WorkerStatus(wx.PyEvent):
    """Simple event that sends text signals to parent for status bar. """
    def __init__(self, string):
        """Init event"""
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_STATUS_ID)
        self.string = string

class WorkerDone(wx.PyEvent):
    """Simple event that signals completion of worker thread."""
    def __init__(self):
        """Init event"""
        wx.PyEvent.__init__(self)
        self.SetEventType(EVT_THREAD_DONE)

class WorkerThread(Thread):
    """ Worker thread class so that GUI stays responsive. Requires the calling
    window to identify itself (notify_window) and send experiment parameters
    (params)."""
    def __init__(self, notify_window, params):
        # Initialize new thread
        Thread.__init__(self)
        
        self._notify_window = notify_window
        self.params = params

        # These identify whether to keep running (_running) and whether or not
        # the parent window has closed (in which case one should not try to post
        # any new events.        
        self._running = True
        self._parent_closing = False

        # Automatically start the thread (could also be done from the gui)        
        self.start()
        
    def run(self):
        """ Defines the tasks for the thread to run in the background of the gui """
        
        # Initialize instruments
        # TODO: Write a config file that contains all the addresses
    #       and configuration information for the entire experiment.
        self.lockin = harmInstr.SR830()
        self.lockin.Initialize()
        #self.filt = harmInstr.FD9002()
        
        # Initialize data structure
        data = dict()
        data['params'] = self.params        # parameters
        data['data'] = list()               # list of dictionaries containing 1st, 2nd, and 3rd harmonic data
        data['pvals'] = dict()              # see following comment
        data['perrs'] = dict()

        # Define appropriate shape for pvals and perrs parameters
        # pvals = vector of best fit values from fitfunc: V0, a, b, c
        # perrs = errors of these parameters, estimated from the covariance matrix
        pvals = np.empty((1,2))
        perrs = np.empty((1,2))
         
        # Until _running is set to False, continue collecting data
        while self._running:
            
            # Collect data
            data['data'].append(self.collectCurve())

            # Calculate slope and offset from X,Y data
            pval, perr, pvalX, pvalY = analyzeXYData(data['data'][-1]['V'],
                            data['data'][-1]['X'], data['data'][-1]['Xerr'],
                            data['data'][-1]['Y'], data['data'][-1]['Yerr'])
                                    
            # Append fit values to pvals array
            pvals = np.append(pvals, [pval], axis=0)
            perrs = np.append(perrs, [perr], axis=0)
            
            # Write values into data structure. This seems a bit slow, but
            # for some reason trying to append values directly to an initialized
            # list adds two values instead of 1            
            data['pvals']['b'] = pvals[1:,0]
            data['pvalX'] = pvalX
#            data['pvals']['a'] = pvals[1:,1]
            data['pvals']['V0'] = pvals[1:,1]
            data['pvalY'] = pvalY
#            data['pvals']['c'] = pvals[1:,3]
            
            data['perrs']['b'] = perrs[1:,0]
#            data['perrs']['a'] = perrs[1:,1]
            data['perrs']['V0'] = perrs[1:,1]
#            data['perrs']['c'] = perrs[1:,3]
            
            # Post data to gui if the gui has not sent a signal that it is closing
            self.PostEvent(ResultEvent(data))
        
        # close instruments when finished
        self.lockin.close()
        #self.filt.close()
        
        # signal completion to parent
        self.PostEvent(WorkerStatus(""))
        self.PostEvent(WorkerDone())
            
    def PostEvent(self, event):
        """ Used to check if parent window exists and if so, post the event. """
        if not self._parent_closing:
            wx.PostEvent(self._notify_window, event)

    def collectCurve(self):
        """ Collect a single harmonic series experiment including capacitance, dc
            voltage sweep for second harmonics, and third harmonic amplitude (single
            amplitude only)
            
            Updated 3/10/2015 to calculate sleep times based on filter time constant
            
            Modified 5/11/2015 for current-only measurement """
        
        params = self.params
        lockin = self.lockin
        
        # set lockin filter settings
        lockin.filter = params['lockinFilter']
        lockin.filterslope = params['lockinFilterSlope']
        lockin.ctrl.write("ISRC 2") # set to perform current measurement with I*10^6 sensitivity
        lockin.sensitivity = params['lockinCapSensitivity']
        lockin.amplitude = params['acAmplitude']/params['acInputGain']

        # get time constant of filter
        tc = lockin.filterdict[lockin.filter][1]
        
        # Set a global scaling factor for sleep times. This comes about because
        # the sleep times were historically set for a 30 ms time constant and should
        # be proportional to the time constant.
        tfactor = 1/0.03
    
        Vs = np.arange(params['MinVoltage'], params['MaxVoltage'] + params['StepVoltage'], params['StepVoltage'])
        data = dict(time=time.time()-params['startTime'], V=Vs, X=np.empty(Vs.size), Xerr=np.empty(Vs.size), Y=np.empty(Vs.size), Yerr=np.empty(Vs.size))
        
        # set frequency
        lockin.frequency = params['acFrequency']
        
        # collect first harmonic (in real units: convert from [V]->[pF])
        self.PostEvent(WorkerStatus("Measuring capacitance..."))
        lockin.harmonic = 1
        lockin.phase = 0
        lockin.reserve = params['lockinReserve']
        #lockin.autoReserve()
        
        time.sleep(1.5*tfactor*tc)
        Cdata = lockin.collectPointsRTheta(5,3*tc)
        data['C'] = Cdata[0]*1e3*1e12/(2*np.pi*lockin.amplitude*params['acInputGain']*lockin.frequency)
        data['Cerr'] = Cdata[2]*1e3*1e12/(2*np.pi*lockin.amplitude*params['acInputGain']*lockin.frequency)
        data['CTheta'] = Cdata[1]
        data['CThetaerr'] = Cdata[3]
        #lockin.phase = Cdata[1]+90
        
        # collect second harmonic
        self.PostEvent(WorkerStatus("Measuring second harmonic amplitude..."))
        lockin.sensitivity = params['lockinSensitivity']
        #if lockin.hasOverload():
        #    lockin.ctrl.write("ISRC 2") # if overloading input at highest sensitivity, back down to lower sensitivity
        lockin.harmonic = 2
        lockin.dcVoltage = Vs[0]/params['dcInputGain']
        lockin.reserve = params['lockinReserve']
        time.sleep(1.5*tfactor*tc)
        #lockin.autoPhase()
        #time.sleep(1.5*tfactor*tc)        
        
        for i,V in enumerate(Vs):
            self.PostEvent(WorkerStatus("Measuring second harmonic amplitude... Vdc = %0.1f mV..." % V))
            lockin.dcVoltage = V/params['dcInputGain']
            time.sleep(0.5*tfactor*tc)
            data['X'][i], data['Y'][i], data['Xerr'][i], data['Yerr'][i] = lockin.collectPointsXY(5, 3*tc)
            
            # convert to real units [A]->[pA]
            data['X'][i] = data['X'][i]*1e12
            data['Xerr'][i] = data['Xerr'][i]*1e12
            data['Y'][i] = data['Y'][i]*1e12
            data['Yerr'][i] = data['Yerr'][i]*1e12

    
        lockin.dcVoltage = 0
        lockin.phase = 0
                        
        # collect third harmonic (in real units: convert [A]->[pA])
        self.PostEvent(WorkerStatus("Measuring third harmonic amplitude..."))
        lockin.harmonic = 3
        time.sleep(1.5*tfactor*tc)
        harm3data = lockin.collectPointsXY(5, 3*tc)
        data['harm3X'] = harm3data[0]*1e12
        data['harm3Y'] = harm3data[1]*1e12
        data['harm3Xerr'] = harm3data[2]*1e12
        data['harm3Yerr'] = harm3data[3]*1e12
        
        return data

    def abort(self):
        """ Receives signal from gui to stop collecting data """
        self._running = False
        
    def flagClose(self):
        """ Receives signal from gui that the gui is closing """
        self._parent_closing = True

""" MAIN PROGRAM SECTION """        

class mainFrame(wx.Frame):
    """ Main control window """
    
    def __init__(self, parent, title):
        """ Initialize window appearance """
        wx.Frame.__init__(self, parent, title=title, size=(1500,1000))

        # Bind to event handler for closing main window
        self.Bind(wx.EVT_CLOSE, self.evtClose)
        
        # Three panels: 0. General status, 1. Measurement status, 2. X,Y cursor position
        self.sb = self.CreateStatusBar(3)
        self.SetStatusWidths([180, -1, -1])
        
        panel = wx.Panel(self)
        
        # Create main gridbagsizer object for entire panel
        sizer = wx.GridBagSizer(4,4)

        # Create input parameter boxes
        boxParams = wx.StaticBox(panel, label="Experiment parameters")
        vbox = wx.StaticBoxSizer(boxParams, wx.VERTICAL)
        
        lblMinVoltage = wx.StaticText(panel, label="Minimum dc voltage (mV)")
        lblMaxVoltage = wx.StaticText(panel, label="Maximum dc voltage (mV)")
        lblStepVoltage = wx.StaticText(panel, label="Step dc voltage (mV)")
        lbldcInputGain = wx.StaticText(panel, label="dc voltage input gain (mV/V)")
        lblAmplitude = wx.StaticText(panel, label="ac rms Amplitude (mV)")
        lblFrequency = wx.StaticText(panel, label="ac Frequency (Hz)")
        lblacInputGain = wx.StaticText(panel, label="ac voltage input gain (mV/V)")
        lbllockinFilter = wx.StaticText(panel, label="lockin filter time constant")
        lbllockinFilterSlope = wx.StaticText(panel, label="lockin filter slope")
        lbllockinReserve = wx.StaticText(panel, label="lockin reserve")        
        lbllockinCapSensitivity = wx.StaticText(panel, label="lockin sensitivity for capacitance measurement")
        lbllockinSensitivity = wx.StaticText(panel, label="lockin sensitivity for harmonics measurement")
                
        self.txtMinVoltage = wx.TextCtrl(panel, name="MinVoltage")
        self.txtMaxVoltage = wx.TextCtrl(panel, name="MaxVoltage")
        self.txtStepVoltage = wx.TextCtrl(panel, name="StepVoltage")
        self.txtdcInputGain = wx.TextCtrl(panel, name="dcInputGain")
        self.txtAmplitude = wx.TextCtrl(panel, name="acAmplitude")
        self.txtFrequency = wx.TextCtrl(panel, name="acFrequency")
        self.txtacInputGain = wx.TextCtrl(panel, name="acInputGain")

        filterkeys = [k for v,k in sorted([(v[0], k) for k, v in harmInstr.filterdict.items()])]
        filterslopekeys = [k for v,k in sorted([(v, k) for k, v in harmInstr.filterslopedict.items()])]
        sensitivitykeys = [k for v,k in sorted([(v, k) for k, v in harmInstr.sensitivitydict.items()])]
        reservekeys = [k for v,k in sorted([(v, k) for k, v in harmInstr.reservedict.items()])]
        reservekeys.append('auto')
        
        self.cmblockinFilter = wx.ComboBox(panel, name="lockinFilter", choices=filterkeys)
        self.cmblockinFilterSlope = wx.ComboBox(panel, name="lockinFilterSlope", choices=filterslopekeys)
        self.cmblockinReserve = wx.ComboBox(panel, name="lockinReserve", choices=reservekeys)        
        self.cmblockinCapSensitivity = wx.ComboBox(panel, name="lockinCapSensitivity", choices=sensitivitykeys)
        self.cmblockinSensitivity = wx.ComboBox(panel, name="lockinSensitivity", choices=sensitivitykeys)        
        
        self.btnLoadParams = wx.Button(panel, label="Load Parameters")
        self.btnSaveParams = wx.Button(panel, label="Save Parameters")
        self.btnSetParams = wx.Button(panel, label="Set Parameters on Lockin")
        
        # Define a list of textboxes for easily reading out parameters
        self.txtboxes = [self.txtMinVoltage, self.txtMaxVoltage, self.txtStepVoltage, self.txtdcInputGain,
                         self.txtAmplitude, self.txtFrequency, self.txtacInputGain,
                         self.cmblockinFilter, self.cmblockinFilterSlope, self.cmblockinCapSensitivity, self.cmblockinSensitivity, self.cmblockinReserve]
        
        vbox.AddMany([(lblMinVoltage), (self.txtMinVoltage, 1, wx.EXPAND),
                    (lblMaxVoltage), (self.txtMaxVoltage, 1, wx.EXPAND),
                    (lblStepVoltage), (self.txtStepVoltage, 1, wx.EXPAND),
                    (lbldcInputGain), (self.txtdcInputGain, 1, wx.EXPAND),
                    (lblAmplitude), (self.txtAmplitude, 1, wx.EXPAND),
                    (lblFrequency), (self.txtFrequency, 1, wx.EXPAND),
                    (lblacInputGain), (self.txtacInputGain, 1, wx.EXPAND),
                    (lbllockinFilter), (self.cmblockinFilter, 1, wx.EXPAND),
                    (lbllockinFilterSlope), (self.cmblockinFilterSlope, 1, wx.EXPAND),
                    (lbllockinReserve), (self.cmblockinReserve, 1, wx.EXPAND),                    
                    (lbllockinCapSensitivity), (self.cmblockinCapSensitivity, 1, wx.EXPAND),                   
                    (lbllockinSensitivity), (self.cmblockinSensitivity, 1, wx.EXPAND)])
                    
        vbox.Add(self.btnLoadParams, proportion=2, flag=wx.EXPAND|wx.TOP, border=10)
        vbox.Add(self.btnSaveParams, proportion=2, flag=wx.EXPAND|wx.TOP, border=10)
        vbox.Add(self.btnSetParams, proportion=2, flag=wx.EXPAND|wx.TOP, border=10)
        
        self.btnLoadParams.Bind(wx.EVT_BUTTON, self.LoadParams)
        self.btnSaveParams.Bind(wx.EVT_BUTTON, self.SaveParams)
        self.btnSetParams.Bind(wx.EVT_BUTTON, self.SetParams)
        
        sizer.Add(vbox, pos=(1,0))
        
        # Add tag controls
        
        boxTags = wx.StaticBox(panel, label="Tags")
        sizerTags = wx.StaticBoxSizer(boxTags, wx.VERTICAL)
        
        self.lstTags = wx.ListCtrl(panel, style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
        self.lstTags.InsertColumn(0, "Time (s)", width=60)
        self.lstTags.InsertColumn(1, "Tag", width=120)
        self.btnAddTag = wx.Button(panel, label="Add new tag")
        self.btnAddTag.Bind(wx.EVT_BUTTON, self.evtAddTag)
        self.btnEditTag = wx.Button(panel, label="Edit selected tag")
        self.btnEditTag.Bind(wx.EVT_BUTTON, self.evtEditTag)
        self.btnDeleteTag = wx.Button(panel, label="Delete selected tag")
        self.btnDeleteTag.Bind(wx.EVT_BUTTON, self.evtDeleteTag)
        
        sizerTags.Add(self.lstTags, 5, wx.EXPAND)
        sizerTags.Add(self.btnAddTag, 1, wx.EXPAND, border=5)
        sizerTags.Add(self.btnEditTag, 1, wx.EXPAND, border=5)
        sizerTags.Add(self.btnDeleteTag, 1, wx.EXPAND, border=5)
        
        sizer.Add(sizerTags, pos=(2,0), span=(1,1), flag=wx.EXPAND|wx.TOP|wx.BOTTOM)
        
        # Add main action buttons
        
        btnbox = wx.BoxSizer(wx.VERTICAL)
        
        self.btnStart = wx.Button(panel, label="Start")
        self.btnStart.Bind(wx.EVT_BUTTON, self.evtPushStart)
        
        self.btnStop = wx.Button(panel, label="Stop")
        self.btnStop.Bind(wx.EVT_BUTTON, self.evtPushStop)
        
        self.btnQuit = wx.Button(panel, label="Quit")
        self.btnQuit.Bind(wx.EVT_BUTTON, self.evtPushQuit)
        
        self.btnZap = wx.Button(panel, label="ZAP")
        self.btnZap.Bind(wx.EVT_BUTTON, self.evtPushZap)
        self.btnZap.SetBackgroundColor = "RED"
        self.btnZap.SetForegroundColor = "WHITE"
        
        btnbox.Add(self.btnStart, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.BOTTOM, border=5)
        btnbox.Add(self.btnStop, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.BOTTOM, border=5)
        btnbox.Add(self.btnQuit, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.BOTTOM, border=5)
        btnbox.Add(self.btnZap, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.BOTTOM, border=5)
        
        sizer.Add(btnbox, pos=(3,0), span=(1,1), flag=wx.EXPAND)
        
        # Add the plots in a GridSizer
        
        gs = wx.GridSizer(2, 3, 10, 10)
                
        self.plt00 = wxmplot.PlotPanel(panel, size=(1.,1.), messenger=self.UpdateStatusBar)
        self.plt01 = wxmplot.PlotPanel(panel, size=(1.,1.), messenger=self.UpdateStatusBar)
        self.plt02 = wxmplot.PlotPanel(panel, size=(1.,1.), messenger=self.UpdateStatusBar)
        self.plt10 = wxmplot.PlotPanel(panel, size=(1.,1.), messenger=self.UpdateStatusBar)
        self.plt11 = wxmplot.PlotPanel(panel, size=(1.,1.), messenger=self.UpdateStatusBar)
        self.plt12 = wxmplot.PlotPanel(panel, size=(1.,1.), messenger=self.UpdateStatusBar)
        
        gs.Add(self.plt00, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT)
        gs.Add(self.plt01, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT)
        gs.Add(self.plt02, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT)        
        gs.Add(self.plt10, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT)
        gs.Add(self.plt11, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT)                
        gs.Add(self.plt12, 1, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT)
        
        sizer.Add(gs, pos=(1,1), span=(2,1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT)
        
        # Initialize the plots. For some reason plt00 has to be plotted twice
        # for labelfontsize to take effect
        self.plt00.plot([0.], [0.], dy=[0.], xlabel = "dc Voltage (mV)", ylabel="2nd harmonic in-phase amplitude (pA)", linewidth=0, marker="o", labelfontsize=6)
        self.plt00.plot([0.], [0.], dy=[0.], xlabel = "dc Voltage (mV)", ylabel="2nd harmonic in-phase amplitude (pA)", linewidth=0, marker="o", labelfontsize=6)
        self.plt01.plot([0.], [0.], dy =[0.], marker='o', xlabel = "dc Voltage (mV)", ylabel="2nd harmonic out-of-phase amplitude (pA)", linewidth=0, labelfontsize=6)
        self.plt02.plot([0], [0], linewidth=0.5, xlabel = "Time (s)", ylabel = "Capacitance (pF)", labelfontsize=6)
        self.plt10.plot([0], [0], linewidth=0.5, xlabel = "Time (s)", ylabel = "2nd harmonic offset (mV)", labelfontsize=6)
        self.plt11.plot([0], [0], linewidth = 0.5, xlabel = "Time (s)", ylabel = "2nd harmonic slope (pA/mV) * 1000", labelfontsize=6)
        self.plt12.plot([0], [0], linewidth = 0.5, xlabel = "Time (s)", ylabel = "3rd harmonic amplitude (pA)", labelfontsize=6)
        
        # Add a status textbox
#        txtStdOut = wx.TextCtrl(panel, style=wx.HSCROLL|wx.VSCROLL|wx.TE_MULTILINE|wx.TE_READONLY)
#        txtStdOut.SetFont(wx.Font(8, wx.MODERN, wx.NORMAL, wx.NORMAL))
#        sizer.Add(txtStdOut, pos=(3,1), span=(1,1), flag=wx.EXPAND|wx.TOP|wx.BOTTOM|wx.RIGHT, border=10)
        
        # Save original stdout references
        self.stdoutref = sys.stdout
        self.stderrref = sys.stderr
        
        # Redirect stdout and stderr to status textbox
#        redirout = RedirectText(txtStdOut)
#        sys.stdout = redirout
#        sys.stderr = redirout

        # Add a data list control
        self.lstData = wx.ListCtrl(panel, style=wx.LC_REPORT|wx.LC_HRULES|wx.LC_VRULES)
        colwidth = 190
        self.lstData.InsertColumn(0, "N", width=40)
        self.lstData.InsertColumn(1, "Time (s)", width=80)
        self.lstData.InsertColumn(2, "C (pF)", width=80)
        self.lstData.InsertColumn(3, "Second harmonic offset (mV)", width=colwidth)
        self.lstData.InsertColumn(5, "Second harmonic slope (pA/mV)", width=colwidth)
        self.lstData.InsertColumn(7, "Third harmonic amplitude (pA)", width=colwidth)
        sizer.Add(self.lstData, pos=(3,1), span=(1,1), flag=wx.EXPAND, border=5)
        
        # Add the file name label to the top of the screen
        self.lblFileName = wx.StaticText(panel, label=" ", style=wx.ALIGN_CENTER)
        self.lblFileName.SetFont(wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.BOLD))
        sizer.Add(self.lblFileName, pos=(0,1), flag=wx.EXPAND)
        
        # Allow things to stretch if the window size is changed
        sizer.AddGrowableCol(1)
        sizer.AddGrowableRow(1)
        sizer.AddGrowableRow(2)
        sizer.AddGrowableRow(3)
        
        panel.SetSizer(sizer)

        # Load in default parameters
        self.LoadParams(None, filename="defaultI.json")
        
        # Set up event handlers for worker thread
        EVT_RESULT(self,self.evtResult)
        EVT_STATUS(self, self.evtWorkerStatus)
        EVT_DONE(self, self.evtWorkerDone)
        self.worker = None
        
        self.sizer = sizer          # So it can be called later
        self.sizer.Layout()
        
        # Enable parameter box, disable tags and stop button
        self.TagsEnabled(False)
        self.EnableParams()
        self.btnStop.Enable(False)
        
        # Show window
        self.Center()
        self.Show(True)
        
    def UpdateStatusBar(self, string, panel=2):
        """Function for updating status bar. Can't directly use SetStatusText
        because of 'panel' argument"""
        self.SetStatusText(string, number=2)
        
    def ReadParams(self):
        """Read the parameters from the parameter textboxes. Note that the keys
           for the parameter dictionary are stored in the Name textbox.
           TODO: create a config file with the parameter names that is read so
                 this can be dynamically changed"""

        params=dict()
        for txtbox in self.txtboxes:
            if txtbox.Name[:6] == 'lockin':
                params[txtbox.Name] = txtbox.GetValue()
            else:
                params[txtbox.Name] = float(txtbox.GetValue())
                    
        return params
    
    def DisableParams(self):
        """ Formats appearance of disabled textboxes and disables buttons """
        for txtbox in self.txtboxes:
            txtbox.SetEditable(False)
            txtbox.SetBackgroundColour((220,220,220))
            
        self.btnLoadParams.Enable(False)
        self.btnSaveParams.Enable(False)
            
        self.Refresh()
        
    def EnableParams(self):
        """ Formats appearance of enabled textboxes and enables buttons """
        for txtbox in self.txtboxes:
            txtbox.SetEditable(True)
            txtbox.SetBackgroundColour("WHITE")
            
        self.btnLoadParams.Enable(True)
        self.btnSaveParams.Enable(True)

        self.Refresh()
        
    def evtAddTag(self, event):
        """ Asks user for new tag text and adds the tag to the end of the list.
            Note that because the time is stored, it doesn't matter if the tags
            get out of order."""
        dialog = wx.TextEntryDialog(self, "Enter new tag text:", caption="")
        
        if (dialog.ShowModal() == wx.ID_OK) & (len(dialog.GetValue()) != 0):
            tagtext = dialog.GetValue()
            idx = self.lstTags.GetItemCount()
            self.lstTags.InsertStringItem(idx, "%.1f" % (time.time()-self.starttime))
            self.lstTags.SetStringItem(idx, 1, tagtext)
            self.lstTags.EnsureVisible(idx)
    
    def evtEditTag(self, event):
        """ User might wish to edit the selected tag. If no tag is selected,
            nothing happens."""
        
        selection = self.lstTags.GetFirstSelected()
        
        if selection != -1:
            dialog = wx.TextEntryDialog(self, "Edit tag text:", defaultValue=self.lstTags.GetItem(selection, 1).GetText())
            
            if (dialog.ShowModal() == wx.ID_OK):
                tagtext = dialog.GetValue()
                self.lstTags.SetStringItem(selection, 1, tagtext)
                            
    def evtDeleteTag(self, event):
        """ User might wish to delete the selected tag. If no tag is selected,
            nothing happens. """
        
        selection = self.lstTags.GetFirstSelected()
        
        if selection != -1:
            
            self.lstTags.DeleteItem(selection)        
            
    def GetTags(self):
        """ Returns a list of dictionaries containing all the tag information.
            Used for saving data and writing text file. """
        
        tags = list()
        for j in range(self.lstTags.GetItemCount()):
            tags.append(dict())
            tags[j]['time'] = self.lstTags.GetItem(j, 0).GetText()
            tags[j]['tag'] = self.lstTags.GetItem(j, 1).GetText()
        
        return tags
    
    def TagsEnabled(self, enabled):
        """ Enables or disables tag buttons """
        
        self.btnAddTag.Enable(enabled)
        self.btnEditTag.Enable(enabled)
        self.btnDeleteTag.Enable(enabled)
    
    def evtPushStart(self, event):
        """ Event handler for "Start" button """
        self.DisableParams()
        self.btnStart.Enable(False)
        self.btnZap.Enable(False)
        self.btnStop.Enable(True)
        self.TagsEnabled(True)
        
        self.SetStatusText("Collecting data...", number=0)
        
        # if starting new experiment, delete all existing tags and data
        self.lstTags.DeleteAllItems()
        self.lstData.DeleteAllItems()

        # launch scanning subroutine (separate function not necessary but cleaner)        
        self.InitializeScan()
        
    def evtPushStop(self, event):
        """ Event handler for "Stop" button """
                
        """ Now handled by evtWorkerDone. """
        #self.EnableParams()
        #self.btnStart.Enable(True)
        #self.TagsEnabled(False)
        # do a screen capture
        #self.ScreenCapture()
        
        self.btnStop.Enable(False)
        
        self.SetStatusText("Terminating data collection...", number=0)
        
        # send stop signal to the worker thread
        self.worker.abort()
    
    def evtPushQuit(self, event):
        """ Event handler for "Quit" button. Note that the event handler 
            evtClose does the same things whether or not Close() is called from
            here or from the user closing the window directly """
        self.Close()
        
    def evtPushZap(self, event):
        """ Event handler for "Zap" button. Zaps the membrane with 1 V for a
            short time. """
            
        lockin = harmInstr.SR830()
        initV = lockin.dcVoltage
        lockin.dcVoltage = 2
        time.sleep(0.01)
        lockin.dcVoltage = initV
        lockin.close()
    
    def evtClose(self, event):
        """ Event handler for the main window closing """
        # Tell the worker thread to stop, and inform it that the main window has closed
        if type(self.worker) != type(None):
            self.worker.flagClose()
            self.worker.abort()
        
        # Redirect stdout and stderr to their original values
        sys.stdout = self.stdoutref
        sys.stderr = self.stderrref
        
        # If a file has been opened (the "Start" button has ever been pushed),
        # and the stop button was NOT pushed, do a screen capture.
        if hasattr(self, "fn") & self.btnStop.IsEnabled():
            self.ScreenCapture()
        
        # Close window. Can't use self.Close() or evtClose is triggered again
        self.Destroy()
        
    def ScreenCapture(self):
        """ Takes a screenshot of the entire window """
        context = wx.ClientDC( self )
        memory = wx.MemoryDC( )
        x, y = self.ClientSize
        bitmap = wx.EmptyBitmap( x, y, -1 )
        memory.SelectObject( bitmap )
        memory.Blit( 0, 0, x, y, context, 0, 0)
        memory.SelectObject( wx.NullBitmap)
        bitmap.SaveFile(self.fn[:-3] + 'png', wx.BITMAP_TYPE_PNG )
    
    def evtResult(self, event):
        """ Event handler for posting of data from worker thread:
                1. Saves a text file with the capacitance, second harmonic slope
                    and offset parameters, and the third harmonic amplitude in 
                    columnar format. Tags and experiment parameters are written
                    as a header.
                2. Saves a json file with all the raw data. Useful if any 
                    re-analysis (e.g. fit to a different fitfunc) is required.
                    Contains all the tags and experiment parameters.
                3. Updates gui plots and status box
        """
        data = event.data
        data['tags'] = self.GetTags()
        
        #print data['V0'], data['slope0']
        
        """ Text file output """
        # Create columnar data array
        exportdata = np.vstack(([data['data'][j]['time'] for j in range(len(data['data']))],
                                [data['data'][j]['C'] for j in range(len(data['data']))],
                                data['pvals']['V0'], np.array(data['pvals']['b']),
                                [data['data'][j]['harm3X'] for j in range(len(data['data']))]
                                ))

        # Create header from parameters, starting with the filename
        # and starting time and date
        header = self.fn[:-3] + "txt\n" + data['params']['startString']+ "\n"
        for key in data['params'].keys():
            if key != "startString":
                header = header + "%s: %s\n" % (key, data['params'][key])

        # Add tags to header
        header = header + "\nTime (s)\tTag\n"
        for j in range(len(data['tags'])):
            header = header + "%s\t%s\n" % (data['tags'][j]['time'], data['tags'][j]['tag'])

        # Add column headings
        header = header + "\nTime (s)\tC (pF)    \tOffset (mV)\tSlope (pA/mV)\t3rd harmonic amplitude (pA)"
        
        # Write columnar data
        np.savetxt(open(self.fn[:-3] + "txt", "w"), np.transpose(exportdata), header=header, delimiter="\t", fmt="%0.6e")
        
        """ JSON file output. Can fail if file is being read, so add error handler """
        try:
            json.dump(data, open(self.fn, "w"), cls=NumpyAwareJSONEncoder)
        except:
            pass

        """ Update plots """
        # First plot: Last second harmonic measurement
        self.UpdatePlot(self.plt00, data['data'][-1]['V'], data['data'][-1]['X'], dy=data['data'][-1]['Xerr'],
        xlabel = "dc Voltage (mV)", ylabel="2nd harmonic in-phase (pA)", linewidth=0, marker="o", labelfontsize=6)
        x = np.arange(data['params']['MinVoltage'], data['params']['MaxVoltage'], (data['params']['MaxVoltage']-data['params']['MinVoltage'])/100)
        self.plt00.oplot(x, fitfunc(x, data['pvalX'][0], data['pvalX'][1]))
        
        # Second plot: Phase
        #self.UpdatePlot(self.plt01, data['data'][-1]['V'], data['data'][-1]['Theta'], dy=data['data'][-1]['Thetaerr'], marker='o',
        #                xlabel = "dc Voltage (mV)", ylabel="2nd harmonic phase (deg)", linewidth=0, labelfontsize=6)
        self.UpdatePlot(self.plt01, data['data'][-1]['V'], data['data'][-1]['Y'], dy=data['data'][-1]['Yerr'],
        xlabel = "dc Voltage (mV)", ylabel="2nd harmonic out-of-phase (pA)", linewidth=0, marker="o", labelfontsize=6)
        self.plt01.oplot(x, fitfunc(x, data['pvalY'][0], data['pvalY'][1]))

        # Third plot: Capacitance
        self.UpdatePlot(self.plt02, [data['data'][j]['time'] for j in range(len(data['data']))],
                        [data['data'][j]['C'] for j in range(len(data['data']))], linewidth=0.5,
                        xlabel = "Time (s)", ylabel = "Capacitance (pF)", labelfontsize=6)
                        
        # Fourth plot: 2nd harmonic offset
        self.UpdatePlot(self.plt10, [data['data'][j]['time'] for j in range(len(data['data']))],
                        data['pvals']['V0'], linewidth=0.5,
                        xlabel = "Time (s)", ylabel = "2nd harmonic offset (mV)", labelfontsize=6)
                        
        # Fifth plot: 2nd harmonic slope.
        self.UpdatePlot(self.plt11, [data['data'][j]['time'] for j in range(len(data['data']))],
                        data['pvals']['b']*1000, linewidth = 0.5,
                        xlabel = "Time (s)", ylabel = "2nd harmonic slope (pA/mV) * 1000", labelfontsize=6)
        
        # Sixth plot: 3rd harmonic amplitude, in-phase component only
        self.UpdatePlot(self.plt12, [data['data'][j]['time'] for j in range(len(data['data']))],
                        [data['data'][j]['harm3X'] for j in range(len(data['data']))], linewidth = 0.5,
                        xlabel = "Time (s)", ylabel = "3rd harmonic in-phase amplitude (pA)", labelfontsize=6)

        # Print to stdout. #DONE: make this prettier, or turn it into a listctrl
#        print len(data['data'])-1, round(data['data'][-1]['time']), data['data'][-1]['C'], [data['pvals']['V0'][-1], data['pvals']['a'][-1], data['pvals']['b'][-1], data['pvals']['c'][-1]], data['data'][-1]['harm3']
        #sys.stdout.flush()    
        
        # Populate lstData control
        idx = self.lstData.GetItemCount()
        self.lstData.InsertStringItem(idx, "%i" % (len(data['data'])-1))
        self.lstData.SetStringItem(idx, 1, "%i" % round(data['data'][-1]['time']))
        self.lstData.SetStringItem(idx, 2, "%0.6g %s %0.6g" % (data['data'][-1]['C'], chr(177), data['data'][-1]['Cerr']))
        self.lstData.SetStringItem(idx, 3, "%0.6g %s %0.6g" % (data['pvals']['V0'][-1], chr(177), data['perrs']['V0'][-1]))
        self.lstData.SetStringItem(idx, 4, "%0.6g %s %0.6g" % (data['pvals']['b'][-1], chr(177), data['perrs']['b'][-1]))
        self.lstData.SetStringItem(idx, 5, "%0.6g %s %0.6g" % (data['data'][-1]['harm3X'], chr(177), data['data'][-1]['harm3Xerr']))
        self.lstData.EnsureVisible(idx)
        self.ResizelstData()
        
    def ResizelstData(self):
        """ Function for resizing lstData control """
        widths = [40, 80, 0, 0, 0, 0]
        totalwidth = self.lstData.GetSize().GetWidth()
        
        for col in range(self.lstData.GetColumnCount()):
            if widths[col] > 0:
                self.lstData.SetColumnWidth(col, widths[col])
            else:
                self.lstData.SetColumnWidth(col, round((totalwidth-sum(widths)-25)/widths.count(0)))
    
    def evtWorkerStatus(self, event):
        """ Event handler for worker status update. """
        self.SetStatusText(event.string, number=1)    
                
    def evtWorkerDone(self, event):
        """ Event handler for when worker is done. """
        
        self.EnableParams()
        self.btnStart.Enable(True)
        self.btnZap.Enable(True)
        self.TagsEnabled(False)
        
        self.SetStatusText("", number=0)
        
        # do a screen capture
        self.ScreenCapture()
        
    def LoadParams(self, event, filename=None):
        """ Load parameters from JSON file """

        if type(filename) == type(None):        
            dialog = wx.FileDialog(self, wildcard="JSON files (*.json)|*.json")
        
            if dialog.ShowModal() == wx.ID_OK:
                filename = dialog.GetPath()
            else:
                return
        
        # Populate textboxes
        params = json.load(open(filename, "r"))
        for txtbox in self.txtboxes:
            txtbox.SetValue(str(params[txtbox.Name]))

    def SaveParams(self, event):
        """ Save parameters to JSON file """
        params = self.ReadParams()
        
        dialog = wx.FileDialog(self, wildcard="JSON files (*.json)|*.json")
        
        if dialog.ShowModal() == wx.ID_OK:
            retfn = dialog.GetPath()
            json.dump(params, open(retfn, "w"), cls=NumpyAwareJSONEncoder)
            
    def SetParams(self, event):
        """ Set parameters on lockin."""
        
        lockin = harmInstr.SR830()
        params = self.ReadParams()
        lockin.frequency = params['acFrequency']
        lockin.amplitude = params['acAmplitude']/params['acInputGain']
        lockin.filter = params['lockinFilter']
        lockin.filterslope = params['lockinFilterSlope']
        lockin.sensitivity = params['lockinSensitivity']
        lockin.reserve = params['lockinReserve']
        lockin.close()
        
    def UpdatePlot(self, plt, x,y, **kwargs):
        """ Clears and re-writes a plot. Passes kwargs through. """
        plt.clear()
        plt.plot(x,y, **kwargs)
        
    def InitializeScan(self):
        """ Function for initializing a data scan. Sets parameters and filenames
            and starts worker thread. """
        
        # Read parameters from text boxes
        harm2params = self.ReadParams()
        self.starttime = time.time()
        harm2params['startTime'] = self.starttime
        #print 'startTime = ' + `harm2params['startTime']`
        
        # Determine appropriate file name        
        fnprefix = time.strftime("%Y%m%d") + "_"
        fnsuffix = ".dat"
        
        # Find last file name with fnprefix
        flist = glob.glob(fnprefix + "???" + fnsuffix)
        if not bool(len(flist)):
            self.fn = fnprefix + "000" + fnsuffix
        else:
            self.fn = fnprefix + "{:03}".format(len(flist)) + fnsuffix
        
        #print "New file: " + self.fn
        startstring = time.strftime("%m/%d/%Y %I:%M:%S %p")
        #print "Started: " + startstring
        harm2params['startString'] = startstring
        #sys.stdout.flush()
        
        # Update GUI with new file name
        self.lblFileName.SetLabel(self.fn + " (" + startstring + ")")
        self.sizer.Layout()
        
        # Start data collection by initializing new worker thread
        self.worker = WorkerThread(self, harm2params)

                        
# Run the app in such a way that file can be used as library
if __name__ == '__main__':
    app = wx.App()
    frame = mainFrame(None, "Second harmonic analysis")
    app.MainLoop()

