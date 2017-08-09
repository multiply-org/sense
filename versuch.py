#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Filename: versuch.py
# Author: "Alexander Löw"
# Date:
# Last Modified by:   "Thomas Weiß"
# Last Modified time: 2017-08-09 11:09:22

"""
aAJHFHASfdjklA
"""

import numpy as np



from sense import model
import matplotlib.pyplot as plt
from sense.soil import Soil
from sense.canopy import OneLayer


theta = np.deg2rad(np.arange(1.,80.))
freq = 5.

stype = 'turbid_rayleigh'
#stype='turbid_isotropic'
models = {'surface': 'Oh92', 'canopy': stype}
eps = 15. -0.j
s=0.02
d = 0.22
ke = 1.
omega=0.2


soil = Soil(eps=eps, f=freq, s=s)
can = OneLayer(ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke)
S = model.SingleScatRT(surface=soil, canopy=can, models=models, theta=theta, freq=freq)
S.sigma0()


f = plt.figure()
ax = f.add_subplot(111)
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0g']['vv']), color='red', linestyle='-', label='vv s0g')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0c']['vv']), color='blue', linestyle='-', label='vv s0c')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0cgt']['vv']), color='green', linestyle='-', label='vv s0cgt')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['stot']['vv']), color='black', linestyle='-', label='vv stot')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0gcg']['vv']), color='green', linestyle='--', label='vv s0gcg')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['stot']['hh']), color='orange', linestyle='--', label='hh stot')
#ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['stot']['hv']), color='purple', linestyle='-', label='hv stot')
# ax.plot(np.rad2deg(self.theta), self.h, color='blue', linestyle='--', label='H')
ax.set_ylim(-14.,-5.)
ax.set_xlim(0.,70.)
ax.grid()
ax.legend(loc = 3)
plt.show()

theta2=np.deg2rad(50)
Sf = model.SingleScatRT(surface=soil, canopy=can, models=models, theta=theta2 , freq=freq)
Sf.sigma0()
print('stot', 10*np.log10(Sf.__dict__['stot']['vv']))
print('s0cgt', 10*np.log10(Sf.__dict__['s0cgt']['vv']))
print('s0c',10*np.log10(Sf.__dict__['s0c']['vv']))
print('s0gcg',10*np.log10(Sf.__dict__['s0gcg']['vv']))
print('s0g', 10*np.log10(Sf.__dict__['s0g']['vv']))


