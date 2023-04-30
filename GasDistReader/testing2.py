import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
import math
import scipy.interpolate as interp
import scipy.optimize as optimize
from scipy.integrate import quad


a1 = 3.30
a2 = 4.35
a3 = 5.32
a4 = 4.75

rmin1 = 2.00
rmin2 = 3.31
rmin3 = 3.89
rmin4 = 3.19

#Constants
thmin1 = 1.05
thmin2 = 2.62
thmin3 = 4.19
thmin4 = 5.76

def ArmAngle(rad):

    # params:


    arm1 = a1*np.log(rad/rmin1) + thmin1
    arm2 = a2*np.log(rad/rmin2) + thmin2
    arm3 = a3*np.log(rad/rmin3) + thmin3
    arm4 = a4*np.log(rad/rmin4) + thmin4

    # My idea to compactify things a bit, then if I want angle of arm i at R I can just do ArmAngle(R)[i]
    return arm1



def MinDist(x, y):

    r = np.sqrt(x**2 + y**2)
    th = math.atan2(y, x)

    def DerivDist(rad):
        r = np.sqrt(x**2 + y**2)
        th = math.atan2(y, x)

        D = th + a1*np.log(rmin1) - thmin1

        return rad - r*(np.cos(D - a1*np.log(rad)) + a1*np.sin(D - a1*np.log(rad)))

    frq = int(r*10)
    g = np.linspace(0.05, r, frq)

    zeros = np.zeros((frq))
    values = np.zeros((frq))
    for i in range(frq):
        zeros[i] = optimize.root_scalar(DerivDist, method="secant", x0=0.049, x1=g[i]).root 
        values[i] = r**2 + zeros[i]**2 - 2*r*zeros[i]*np.cos(th - ArmAngle(zeros[i]))

    mindist = min(values)

    return mindist

numgrid = 100
a = np.linspace(8, 10, numgrid)

z = np.zeros((numgrid))
for i in range(numgrid):
    z[i] = MinDist(0, a[i])

plt.plot(a, z)
plt.show()




    # x = np.linspace(0.1, y, 20)
    # roots = np.zeros((20))
    # for i in range(20):
        # roots[i] = optimize.root_scalar(DerivDist, method="newton", x0 = x[i], fprime=True)


    

"""
a = np.linspace(0, 7, 100)

z = np.zeros((100))
for i in range(100):
    z[i] = MinDist(a[i])[0]


plt.plot(a, z)

plt.show()
"""
"""

# Model for the profile perpendicular to the center of the arms, a scaled Gaussian:

def ArmScale(x, y):

    r = np.sqrt(x**2 + y**2)
    th = math.atan2(y, x)

    # Width of Gaussian
    sg = 0.6

    # Will find minimal distance squared to arm 1 and then figure out the scaling factor:
    def dist1(rad):
        d1 = r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle(rad))
        return d1

    di = optimize.minimize_scalar(Distance)

  
    # These are gonna be the scaling factors
    return (1/(sg*np.sqrt(2*np.pi)))*np.e**(-0.5*(float(np.sqrt(di.fun)/sg)**2))



numgrid = 100

a = np.linspace(-5, 5, numgrid)
b = np.linspace(-5, 5, numgrid)

zz = np.zeros((numgrid, numgrid))
for m in range(numgrid):
    for n in range(numgrid):
        point = (a[n], b[m])
        zz[m][n] = ArmScale(a[n], b[m])
# print(a[47])
# print(b[37])
print(ArmScale(-0.2525252525252526, -1.2626262626262625))
print(zz[37][47])


h = plt.pcolormesh(a, b, zz, shading='auto')

plt.plot(-0.25, -1.2, marker="o", markersize=5, markerfacecolor="white")

plt.axis('scaled')
plt.colorbar()
plt.xlabel("Kpc")
plt.ylabel("Kpc")
plt.title("H1 Surface Density")

plt.show()


"""