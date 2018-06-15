"""
Fig 4.31 in Ulaby (2014)
check of Dobson 85 model

Fig 4.31 are measurement data not model data!

Shape of the graphs looks fine
different bulk densities result in different graphs!
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)) + os.sep + '..')
from sense.dielectric import Dobson85


mv = np.linspace(0.,0.5, 100)
f = 5.

fig = plt.figure()
ax1 = fig.add_subplot(111)

D = Dobson85(sand=0.515, clay=0.135,freq=f,mv=mv)
# Eqs. 4.69 simplistic approach
e = D.eps

ax1.plot(mv, np.real(e), color='blue')
ax1.plot(mv, np.imag(e), color='red')
ax1.grid()
ax1.set_title('Fig 4-31')
ax1.set_xlim(0.,0.45)
ax1.set_ylim(0.,27.)

D = Dobson85(sand=0.515, clay=0.135,freq=f,mv=mv, debye=True)
# Eqs. 4.67 Debye model with conductivity term for e2
e = D.eps

ax1.plot(mv, np.real(e), color='green')
ax1.plot(mv, np.imag(e), color='green')

D = Dobson85(sand=0.515, clay=0.135,freq=f,mv=mv, debye=True, single_debye=True)
# Eqs. 4.14 single Debye dielectric model for pure water
e = D.eps

ax1.plot(mv, np.real(e), color='black')
ax1.plot(mv, np.imag(e), color='black')


plt.show()
