# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""
compare own results against references
from the Ulaby example codes provided
http://mrs.eecs.umich.edu/codes/Module10_5/Module10_5.html

for the Oh92 model (PRISM1)
"""
import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__)) + os.sep + '..')

from sense.surface import Oh92
from sense.surface import Dubois95
from sense.util import f2lam

import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

eps = 15-0.j   # note that normally the imagionary part is supposed to be negative!

freq = 5.  # GH
s = 0.02  # m
ks = (2.*np.pi/f2lam(freq))*s

theta = np.deg2rad(np.arange(0.,71.) )

O = Oh92(eps, ks, theta)
O.plot()

d =  Dubois95(eps, ks, theta, f2lam(freq))
d.plot()

o = Oh92(eps, ks, np.deg2rad(50))
print(10*np.log10(o.vv))
print(10*np.log10(o.hh))
print(10*np.log10(o.hv))

d = Dubois95(eps, ks, np.deg2rad(40), f2lam(freq))
print(d.vv)
print(d.hh)
