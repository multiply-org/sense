"""
Fig 4.1 in Ulaby (2014)
check of pure water single debye dielectric model

looks pretty good!
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)) + os.sep + '..')
from sense.dielectric import Dobson85


mv = 0.1
f = np.linspace(1.,50., 100)
t = 0.

fig = plt.figure()
ax1 = fig.add_subplot(111)

D = Dobson85(sand=0.5, clay=0.33,freq=f,mv=mv, debye=True, single_debye=True, temp=t)
e = D.ew

ax1.plot(f, np.real(e), color='blue')
ax1.plot(f, np.imag(e), color='red')
ax1.grid()
ax1.set_title('Fig 4-1b')
ax1.set_xlim(1.,50.)
ax1.set_ylim(1.,100.)
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.show()
