"""
compare own results against references
from the Ulaby example codes provided
http://mrs.eecs.umich.edu/codes/Module4_7/Module4_7.html

check of Dobson 85 model
fit is really good!

Difference in complexity of Dobson model only for imaginary part
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)) + os.sep + '..')
from sense.dielectric import Dobson85


mv = 0.3
f = np.linspace(1.,20, 100)

fig = plt.figure(figsize=(8,2))

ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)


D = Dobson85(sand=0.5, clay=0.33,freq=f,mv=mv, bulk=1.7)
# Eqs. 4.69 simplistic approach
e = D.eps

ax1.plot(f, np.real(e), color='blue')
ax2.plot(f, np.imag(e), color='red')
ax1.grid()
ax1.set_title('Fig 4-31')
ax1.set_xlim(1.4,20)
ax1.set_ylim(12.,22.)

D = Dobson85(sand=0.5, clay=0.33,freq=f,mv=mv, debye=True, bulk=1.7)
# Eqs. 4.67 Debye model with conductivity term for e2
e = D.eps

ax1.plot(f, np.real(e), color='green')
ax2.plot(f, np.imag(e), color='green')

D = Dobson85(sand=0.5, clay=0.33,freq=f,mv=mv, debye=True, single_debye=True, bulk=1.7)
# Eqs. 4.14 single Debye dielectric model for pure water
e = D.eps

ax1.plot(f, np.real(e), color='black')
ax2.plot(f, np.imag(e), color='black')

ax2.grid()
ax2.set_title('Fig 4-31')
ax2.set_xlim(1.4,20)
ax2.set_ylim(1.,7.)


plt.show()
