import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
import math
import scipy.interpolate as interp
import scipy.optimize as optimize
from scipy.integrate import quad


# Disk and arm normalizations:

ndisk = 0.760479
narm1 = 0.526987
narm2 = 6.10118e-08
narm3 = 3.1493
narm4 = 8.3478


# Exponential disk with central hole
# Used just for CO

def ExpDisk(rad, norm):

    # Filler values for the free parameters at the moment
    rs = 8 # normalization constant
    r0 = 1.12804
    rh = 6.43578
    hi = 7.30415

    return norm * np.e**(-(rad - rs)/r0) * (1 - np.e**(-(rad/rh)**hi))


def WarpSpline(rad, mode):

    # Filler values for radii, three different modes
    R0_0 = -0.0631513
    R5_0 = -0.00767916
    R10_0 = -0.028553
    R15_0 = 0.0562734
    R20_0 = 0.764476
    R50_0 = 20

    R0_1 = -0.117265
    R5_1 = 0.0517391
    R10_1 = -0.101549
    R15_1 = -0.737641
    R20_1 = -1.69437
    R50_1 = -20

    R0_2 = 0.26564
    R5_2 = -0.0213707
    R10_2 = -0.00681474
    R15_2 = -0.00494624
    R20_2 = 0.570304
    R50_2 = 15.5518

    rvals = np.array([0, 5, 10, 15, 20, 50])
    nvals0 = np.array([R0_0, R5_0, R10_0, R15_0, R20_0, R50_0])
    nvals1 = np.array([R0_1, R5_1, R10_1, R15_1, R20_1, R50_1])
    nvals2 = np.array([R0_2, R5_2, R10_2, R15_2, R20_2, R50_2])

    spline0 = interp.CubicSpline(rvals, nvals0)
    spline1 = interp.CubicSpline(rvals, nvals1)
    spline2 = interp.CubicSpline(rvals, nvals2)

    if mode == 0:
        return spline0(rad)
    elif mode == 1: 
        return spline1(rad)
    elif mode == 2:
        return spline2(rad)       
    else:
        return print("modes only 0, 1, or 2")

def DiskWarp(rad, th):

    # Finding the amplitudes
    w0 = WarpSpline(rad, 0)
    w1 = WarpSpline(rad, 1)
    w2 = WarpSpline(rad, 2)

    # Other free params:
    th1 = 1.47548
    th2 = 2.76487

    return w0 + w1*np.sin(th - th1) + w2*np.sin(2*th - th2)

def Zprime(z, rad, th):

    # Finding z0, the distance from the disk's central plane
    z0 = DiskWarp(rad, th)

    return z - z0



def ScaleHeight(rad):

    # Params:
    rz0 = 8.5 # constant
    zs = 9.90944
    rz = 6.93695

    return zs * np.e**((rad - rz0)/rz)

# For CO, square of hyperbolic secant:

def SqrHypSec(x):
    return (1/np.cosh(x))**2


# In the CO best fit model there's also a central bulge/bar component that we need to add

def CentBulge(X, Y, Z):

    # params:
    nb = 145.332
    ei = 0.442983
    pi = 0.711968
    rb = 0.100744
    zb = 0.00116806
    x0 =-0.578175
    thb = -np.pi / 6 # constant, the paper has -30 degrees idk why they're not using radians


    Xprime = X*np.cos(thb) + Y*np.sin(thb) + x0
    Yprime = -X*np.sin(thb) + Y*np.cos(thb)
    Rprime = np.sqrt(Xprime**2 + (Yprime / 0.3)**2)
    Rr = ((Rprime/rb) + (Z/zb))**(-1)

    return nb * np.e**(-Rr**ei) * Rr**pi



def ArmAngle(rad):

     # params:
    a1 = 3.33783
    a2 = 4.36017
    a3 = 5.36672
    a4 = 4.72076

    rmin1 = 2.03805
    rmin2 = 3.31587
    rmin3 = 3.94189
    rmin4 = 3.16184

    #Constants
    thmin1 = 1.05
    thmin2 = 2.62
    thmin3 = 4.19
    thmin4 = 5.76

    arm1 = a1*np.log(rad/rmin1) + thmin1
    arm2 = a2*np.log(rad/rmin2) + thmin2
    arm3 = a3*np.log(rad/rmin3) + thmin3
    arm4 = a4*np.log(rad/rmin4) + thmin4

    # My idea to compactify things a bit, then if I want angle of arm i at R I can just do ArmAngle(R)[i]
    return (0, arm1, arm2, arm3, arm4)


# Model for the profile perpendicular to the center of the arms, a scaled Gaussian:

def ArmScale(r, th, z):

    # Width of Gaussian
    sg = 0.6

    # Will find minimal distance to arm 1 and then figure out the scaling factor:
    def dist1(rad):
        d1 = np.sqrt(r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle(rad)[1]) + z**2)
        return d1
    def dist2(rad):
        d2 = np.sqrt(r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle(rad)[2]) + z**2)
        return d2
    def dist3(rad):
        d3 = np.sqrt(r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle(rad)[3]) + z**2)
        return d3
    def dist4(rad):
        d4 = np.sqrt(r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle(rad)[4]) + z**2)
        return d4

    if r < 9 + 1e-6:
        mins1, mins2, mins3, mins4 = np.zeros((9)), np.zeros((9)), np.zeros((9)), np.zeros((9))
        for i in range(9):
            mins1[i] = optimize.minimize_scalar(dist1, bounds=(1e-6+i, 1e-6+i+1)).fun
            mins2[i] = optimize.minimize_scalar(dist2, bounds=(1e-6+i, 1e-6+i+1)).fun
            mins3[i] = optimize.minimize_scalar(dist3, bounds=(1e-6+i, 1e-6+i+1)).fun
            mins4[i] = optimize.minimize_scalar(dist4, bounds=(1e-6+i, 1e-6+i+1)).fun
        di = (0, min(mins1), min(mins2), min(mins3), min(mins4))

    else:
        di = (0, optimize.minimize_scalar(dist1, bounds=(1e-6, 30)).fun, optimize.minimize_scalar(dist2, bounds=(1e-6, 30)).fun, optimize.minimize_scalar(dist3, bounds=(1e-6, 30)).fun, optimize.minimize_scalar(dist4, bounds=(1e-6, 30)).fun)

  
    def Gauss(dist):
        return (1/(sg*np.sqrt(2*np.pi)))*np.e**(-0.5*(float(np.sqrt(dist)/sg)**2))
    
    # These are gonna be the scaling factors
    return (0, Gauss(di[1]), Gauss(di[2]), Gauss(di[3]), Gauss(di[4]))


# We now want the final disk model in Cartesian coords
# For HI, best fit is the cubic spline radial profile and sqrt of hyp sec for vertical profile

# For CO, best fit is the exp dist for radial profile and the square of hyp sec for the vertical profile

def DiskNumDensityCO(X, Y, Z, norm):

    Rad = np.sqrt(X**2 + Y**2)
    Theta = np.arctan(Y/X)

    return ExpDisk(Rad, norm) * SqrHypSec(Zprime(Z, Rad, Theta) / ScaleHeight(Rad))





def CODensity(X, Y, Z):

    Rad = np.sqrt(X**2 + Y**2)
    Th = math.atan2(Y, X)

    nddisk = DiskNumDensityCO(Rad, Th, Z, ndisk)
    ndbulge = CentBulge(X, Y, Z)
    ndarm1 = ArmScale(Rad, Th, Z)[1] * DiskNumDensityCO(Rad, Th, Z, narm1)
    ndarm2 = ArmScale(Rad, Th, Z)[2] * DiskNumDensityCO(Rad, Th, Z, narm2)
    ndarm3 = ArmScale(Rad, Th, Z)[3] * DiskNumDensityCO(Rad, Th, Z, narm3)
    ndarm4 = ArmScale(Rad, Th, Z)[4] * DiskNumDensityCO(Rad, Th, Z, narm4)

    return nddisk + ndbulge + ndarm1 + ndarm2 + ndarm3 + ndarm4


def PlaneCODensity(x, y):
    return CODensity (x, y, 0)


print(PlaneCODensity(1,1))


numgrid = 100

a = np.linspace(-10, 10, numgrid)
b = np.linspace(-10, 10, numgrid)

zz = np.zeros((numgrid, numgrid))
for m in range(numgrid):
    for n in range(numgrid):
        point = (a[n], b[m])
        zz[m][n] = PlaneCODensity(a[n], b[m])

np.savetxt("CO_Surface_Density.dat", zz)

h = plt.pcolormesh(a, b, zz, shading='auto')

plt.axis('scaled')
plt.colorbar()
plt.xlabel("Kpc")
plt.ylabel("Kpc")
plt.title("CO Surface Density")

plt.show()
