import numpy as np
import matplotlib
matplotlib.rcParams['backend'] = 'wxagg'
import matplotlib.pyplot as plt
import scipy.ndimage as sn

""" Test effect of Gaussian vs. Boxcar filter on mean """

def autocorr(x):
    return np.correlate(x,x)[0]/len(x)

if 0:

    t = np.arange(0,10000)
    x = np.random.randn(len(t))
    #x[5001] = 100
    fwidth = 10
    xf = sn.uniform_filter1d(x, fwidth)
    xfg = sn.gaussian_filter1d(x, fwidth)

    print np.var(x), np.var(x)/np.var(xf), np.var(x)/np.var(xfg)
    print np.sqrt(np.var(x)/np.var(xf)/fwidth), np.sqrt(np.var(x)/np.var(xfg)/fwidth)
    print np.std(x), np.std(xf)*np.sqrt(fwidth), np.std(xfg)*np.sqrt(fwidth)*1.88

    res = np.correlate(xfg,xfg, 'full')/len(x)
    res = res[(res.size/2):]
    print np.sum(res), res[0], np.correlate(xfg,xfg)/len(x)
    
    print np.std(x)/np.sqrt(autocorr(x)), np.std(xf)/np.sqrt(autocorr(xf)), np.std(xfg)/np.sqrt(autocorr(xfg))
    
""" Result:
        Filtering (and subsequent oversampling) causes the average to be estimated
        too precisely; for a boxcar average, the standard deviation of the averaged 
        data should be multiplied by the root of the boxcar width; for a Gaussian with width
        sigma, by the root of sigma x 1.88."""

""" Now determine correct scaling for covariance matrix """

def mb2V0(m,b):

    return -b/m

def mb2V0_jac(m,b):

    return np.array([b/m**2, -1/m])

if 1:

    x = np.arange(-50.0, 50.0, 0.1)
    y = np.polyval([0.242, 3], x) + 1*np.random.randn(len(x))
    y = sn.gaussian_filter1d(y, 20)

    p0, pcov0 = np.polyfit(x,y,1, full=False, cov='unscaled')
    p0sc, pcov0sc = np.polyfit(x,y,1, full=False, cov=True)    
    # pcov0 is already scaled by chi-squared here, so only need to scale by gaussian filter factor
    # we know sampling rate (say 76 Hz), and we know the time constant (say 100 ms -- use autocorrelation or PSD to calculate any pre/postfactors)
    # Thus we know oversampling, and can correct error bar accordingly. Must be some way to do this based on autocorrelation function.
    print pcov0, pcov0sc, np.sqrt(np.diag(pcov0sc))/p0sc

#    plt.plot(x,y, x, np.polyval(p0, x))
#    plt.show()
    yp = y-np.polyval(p0,x)
    print np.std(yp)**2 / 2
    res = np.correlate(yp, yp, 'full')/len(y)
    res = res[(res.size/2):]
    print np.sum(res), autocorr(yp)
    #plt.plot(res)
    #plt.show()

    normfactor = np.std(yp)**2 / 2

    """ Final piece: figure out how to use Jacobian for transformation """
    J = mb2V0_jac(*p0)
    V0_err = np.sqrt(np.dot(J, np.dot(pcov0, J.T)))
    # Error bar for V0

    #show that in this case root sum of squares is equal to relative error on V0
    print np.sqrt(np.sum((np.sqrt(np.diag(pcov0sc))/p0)**2))
    print V0_err/mb2V0(*p0)
