import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
import math
import scipy.interpolate as interp
import scipy.optimize as optimize
from scipy.integrate import quad
"""
def ArmAngle1(rad):

    # params:
    a1 = 3.33783
    rmin1 = 2.03805
    thmin1 = 1.05

    arm1 = a1*np.log(rad/rmin1) + thmin1

    # My idea to compactify things a bit, then if I want angle of arm i at R I can just do ArmAngle(R)[i]
    return arm1


# Distance from point (x,y) to point (rad, ArmAngle(rad)) in the plane z = 0
def dist1(rad, x, y):
    r = np.sqrt(x**2 + y**2)
    th = math.atan2(y, x)
    d1 = r**2 + rad**2 - 2*r*rad*np.cos(th - ArmAngle1(rad))
    
    return d1



a = np.linspace(0, 5, 100)

z = np.zeros((100))
for i in range(100):
    z[i] = dist1(a[i], -0.25, -1.2)


plt.plot(a, z)

plt.show()

for i in range(10):
    print(i)
"""

a, b, c = 1, 2, 3

print(a)
print(b)
print(c)