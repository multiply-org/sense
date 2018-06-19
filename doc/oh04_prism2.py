"""
comparison with figure 3 in
Oh (2004): Quantitative retrieval of soil moisture content and surface roughness from multipolarized radar observations of bare soil surface. IEEE TGRS 42(3). 596-601.

Oh04 implementation seems to be OK!
"""

import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__)) + os.sep + '..')

import numpy as np

from sense.surface import I2EM

import matplotlib.pyplot as plt

from sense.surface import Oh92, Oh04
from sense.util import f2lam

def db(x):
    return 10.*np.log10(x)


plt.close('all')

f = plt.figure()
ax1 = f.add_subplot(121)
ax2 = f.add_subplot(122)

f  = 5.3 # GHz
lam = f2lam(f)  # m
k = 2.*np.pi/lam
mv = np.linspace(0.01,0.4)
theta = np.deg2rad(40)

s = 0.3/100.
ks = k * s
Oh = Oh04(mv, ks, theta)
ax1.plot(db(Oh.hv), mv)
ax2.plot(db(Oh.p), mv)
s = 1.2/100.
ks = k * s
Oh = Oh04(mv, ks, theta)
ax1.plot(db(Oh.hv), mv)
ax2.plot(db(Oh.p), mv)
s = 4.8/100.
ks = k * s
Oh = Oh04(mv, ks, theta)
ax1.plot(db(Oh.hv), mv)
ax2.plot(db(Oh.p), mv)

ax1.grid()
ax1.set_xlim(-34.,-14.)
ax1.set_ylim(0.01,0.4)
ax1.set_xticks(np.arange(-34, -12, 2))
ax2.grid()
ax2.set_xlim(-4.,1.)
ax2.set_ylim(0.01,0.4)
ax2.set_xticks(np.arange(-4, 1, 0.5))
plt.show()

