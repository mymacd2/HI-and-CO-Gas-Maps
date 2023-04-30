import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
import scipy.interpolate as interp
import scipy.optimize as optimize


def Spline(rad):

    # Filler values for set radius parameter values, RX is X kpc from the Solar System
    # The values calculated in the model are normalized to R8 (so e.g. for the arms this value is gonna be different)

    R8 = 0.521
    R8v = np.e**(R8)

    R0 = np.e**(1.72e-6 * R8)
    R2 = np.e**(0.284 * R8)
    R4 = np.e**(0.807 * R8)
    R6 = np.e**(1.19 * R8)
    R10 = np.e**(0.798 * R8)
    R15 = np.e**(0.477 * R8)
    R20 = np.e**(0.0457 * R8)
    R50 = np.e**(1.12e-3 * R8)
    
    # the first log idea, but this gave me even more negative values lol so I don't think it's the solution
    """
    R0 = np.log(1.72e-6 * R8)
    R2 = np.log(0.284 * R8)
    R4 = np.log(0.807 * R8)
    R6 = np.log(.19 * R8)
    R10 = np.log(0.798 * R8)
    R15 = np.log(0.477 * R8)
    R20 = np.log(0.0457 * R8)
    R50 = np.log(1.12e-3 * R8)
    """

    # Interpolating the cubic spline from the points given
    rvals = np.array([0, 2, 4, 6, 8, 10, 15, 20, 50])
    nvals = np.array([R0, R2, R4, R6, R8v, R10, R15, R20, R50])

    spline = interp.CubicSpline(rvals, nvals)

    return spline(rad)


def SqrtHypSec(x):
    return np.sqrt(1/np.cosh(x))


numgrid = 100

a = np.linspace(0, 50, numgrid)

plt.plot(a, Spline(a))

plt.show()