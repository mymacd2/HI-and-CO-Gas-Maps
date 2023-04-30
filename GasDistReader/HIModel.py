import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
import math
import scipy.interpolate as interp
import scipy.optimize as optimize
from scipy.integrate import quad

# Goal: n(origin) = 3.77e-6)

"""
1) Retranscribe with actual xml numbers
2) Use different optimization function, maybe look at larger steps to make sure I'm getting the real absolute minimum?
"""


"""
This is what I have so far for the paper's HI 3D gas model. Main problems/questions at the moment:

1) Broadly, is the code even correctly transcribing the paper? I think that for the most part it is, a lot of it was just copying down
equations, but in a few areas I'm still on the fence, especially with the spline business. The meaning of "cubic spline in logarithm of number density"
is a bit of a mystery to me still. I think I understand the calculation of the arm densities, but I'm also not certain about this.

2) I'm getting negative values for the density at some places, which is definitely not right--I'm not exactly sure where they're showing up,
but I'm pretty confident that it's in the spline, since the NumDensityHI function is just the spline multiplied by hyperbolic secant,
which is always > 0. So the negative numbers have to be coming from the spline, I think.

3) There's a big speed problem. The code is running really, really slowly--like it's taking me 5+ minutes to run the integrated surface
density charts on only 100 points. I'm not sure exactly where the biggest areas of inefficiency are--I do have a strong suspicion that 
the way I'm plotting the data, with the nested for loops, has got to be very inefficient, but I'm not sure of a better way to do it, since 
the density functions don't take arrays and I'm not sure how to isolate points from a mesh grid without looping through them.
I also bet that the arm density calculations, where I have to minimize distance to each arm, take up a lot of time. I don't know how I would
get around this, however--is there a faster way, maybe by direct calculation, of finding the minimum distance from a point to the arm function?
I tried figuring this out about a week ago but got pretty lost. I'm sure it's possible, and I'll keep working on it, but if you guys have any
tips (or ways to make the minimization go faster) then that would also be very helpful.

4) Looking at the graphs of the z = 0 Galactic plane density I have so far, it seems like (apart from the negative numbers) there's some funny business
going on towards the GC. No real idea of where this is coming from, but I'll keep debugging and looking through my code and hopefully spot what's 
happening.
"""



# Normalizations to be used later

R8disk = 0.63291
R8arm1 = 0.745405
R8arm2 = 0.786137
R8arm3 = 1.29524
R8arm4 = 2.09223


# Cubic spline in logarithm of number density (not 100% sure what "in logarithm of number density" means,
# does it just mean that R0 = log(1.72e-6 * R8), etc? Like we're doing the spline over the logs of the density values
# at the distances given?)
# Used just for HI:

def Spline(rad, norm):

    # Filler values for set radius parameter values, RX is X kpc from the Solar System
    
    R0 = np.log(0.00140987)
    R2 = np.log(0.174074)
    R4 = np.log(0.202369)
    R6 = np.log(0.29917)
    R8 = np.log(0.25)
    R10 = np.log(0.200351)
    R15 = np.log(0.120302)
    R20 = np.log(0.0115289)
    R50 = np.log(0.00028641)

    # Interpolating the cubic spline from the points given
    rvals = np.array([0, 2, 4, 6, 8, 10, 15, 20, 50])
    nvals = np.array([R0, R2, R4, R6, R8, R10, R15, R20, R50])

    spline = interp.CubicSpline(rvals, nvals)

    return norm * np.e**(spline(rad))


# Galactic disk warp: uses another spline to model amplitude radial dependence. 
# Mode is either 0, 1, or 2

def WarpSpline(rad, mode):

    # Filler values for radii, three different modes
    
    R0_0 = -0.0631513
    R5_0 = -0.00767916
    R10_0 = -0.0285536
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
    
    """
    R0_0 = -0.0756
    R5_0 = -0.00819
    R10_0 = -0.0288
    R15_0 = 0.0576
    R20_0 = 0.767
    R50_0 = 20

    R0_1 = -0.146
    R5_1 = 0.0520
    R10_1 = -0.101
    R15_1 = -0.737
    R20_1 = -1.71
    R50_1 = -20

    R0_2 = 0.287
    R5_2 = -0.0192
    R10_2 = -0.00716
    R15_2 = -0.00587
    R20_2 = 0.587
    R50_2 = 14.9
    """

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
    # th1 = 1.47548
    # th2 = 2.76487
    th1 = 4.61
    th2 = 2.73

    return w0 + w1*np.sin(th - th1) + w2*np.sin(2*th - th2)

print(DiskWarp(0, 0))


# Vertical profile model: function of distance from central plane of disk divided by
# scale height of the disk

def Zprime(z, rad, th):

    # Finding z0, the distance from the disk's central plane
    z0 = DiskWarp(rad, th)

    return z - z0

# Function defining scale height as a function of radius:
def ScaleHeight(rad):

    # Params:
    rz0 = 8.5 # constant
    zs = 10.6245 # this feels wrong potentially...
    rz = 6.93695

    return zs * np.e**((rad - rz0)/rz)

print(ScaleHeight(0))

# Functions describing vertical profile of the disk, as a function of rad and Zprime (here, x is Zprime/ScaleHeight)

# For HI, square root of hyperbolic secant:

def SqrtHypSec(x):
    return np.sqrt(1/np.cosh(x))


# Model for the arms, just logorithmic:

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

    # Will find minimal distance squared to arm 1 and then figure out the scaling factor:
    def dist1(rad):
        d1 = r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle(rad)[1]) + z**2
        return d1
    def dist2(rad):
        d2 = r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle(rad)[2]) + z**2
        return d2
    def dist3(rad):
        d3 = r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle(rad)[3]) + z**2
        return d3
    def dist4(rad):
        d4 = r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle(rad)[4]) + z**2
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

def NumDensityHI(Rad, Th, Z, norm):
    return Spline(Rad, norm) * SqrtHypSec(Zprime(Z, Rad, Th) / ScaleHeight(Rad))



##########################################################################################################################################################

# The big boy function, should (hopefully) be returning the number density of HI gas at the point (X, Y, Z)

def HIDensity(X, Y, Z):

    Rad = np.sqrt(X**2 + Y**2)
    Th = math.atan2(Y, X)

    nddisk = NumDensityHI(Rad, Th, Z, R8disk)
    ndarm1 = ArmScale(Rad, Th, Z)[1] * NumDensityHI(Rad, Th, Z, R8arm1)
    ndarm2 = ArmScale(Rad, Th, Z)[2] * NumDensityHI(Rad, Th, Z, R8arm2)
    ndarm3 = ArmScale(Rad, Th, Z)[3] * NumDensityHI(Rad, Th, Z, R8arm3)
    ndarm4 = ArmScale(Rad, Th, Z)[4] * NumDensityHI(Rad, Th, Z, R8arm4)

    return nddisk + ndarm1 + ndarm2 + ndarm3 + ndarm4


def DiskNumDensityHI(X,Y):
    Rad = np.sqrt(X**2 + Y**2)
    Th = math.atan2(Y, X)

    nddisk = NumDensityHI(Rad, Th, 0, R8disk)

    return nddisk

def Arm1NumDensity(X, Y):
    Rad = np.sqrt(X**2 + Y**2)
    Th = math.atan2(Y, X)

    ndarm1 = ArmScale(Rad, Th, 0)[1] * NumDensityHI(Rad, Th, 0, R8arm1)

    return ndarm1


# The density at z = 0, I'm using this to create some quick guiding charts in order to spot bugs (although it's 
# still not all that quick to create them, a few mins for numgrid = 100)

def PlaneHIDensity(x, y):
    return HIDensity(x, y, 0)


# The integrated surface density. This runs very, very slowly––numgrid = 10 finished running in around the same time as
# numgrid = 100 did for the PlaneHIDensity function

def SurfaceHIDensity(x,y):
    def vertd(z):
        return HIDensity(x, y, z)
    return quad(vertd, -5, 5)[0]

x = 0
y = 0
print(PlaneHIDensity(x,y))


"""

# Plotting attempts: definitely very inefficient, any help would be much appreciated

numgrid = 100

a = np.linspace(-5, 5, numgrid)
b = np.linspace(-5, 5, numgrid)

zz = np.zeros((numgrid, numgrid))
for m in range(numgrid):
    for n in range(numgrid):
        point = (a[n], b[m])
        zz[m][n] = PlaneHIDensity(a[n], b[m])

np.savetxt("Surface_Density.dat", zz)

h = plt.pcolormesh(a, b, zz, shading='auto')

plt.axis('scaled')
plt.colorbar()
plt.xlabel("Kpc")
plt.ylabel("Kpc")
plt.title("H1 Surface Density")

plt.show()

"""