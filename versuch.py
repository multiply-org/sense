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

from sense import surface

# from sense.surface import prevot1993
from sense.surface import oh1992

from sense.surface.prevot1993 import Prevot93_surface
from sense.vegetation.prevot1993 import Prevot93_vegetation


from sense import model
from sense.model import Ground

import matplotlib.pyplot as plt

import numpy as np
from sense import model
from sense.soil import Soil
from sense.canopy import OneLayer
from sense.surface import Dubois95
from sense.util import f2lam


# xxx=SingleScatRT(theta=25, surface={"eps": 3}, canopy='turbid_isotropic', models={"surface": 'Oh92', "canopy": 'turbid_isotropic'}, freq=3.)

# theta = np.deg2rad(np.arange(5.,80.))
# theta = np.deg2rad(25.)
# freq = 5.



# stype='turbid_isotropic'
# models = {'surface': 'Oh92', 'canopy': stype}
# eps = 5. -3.j
# soil = Soil(eps=eps, f=5., s=0.02)
# can = OneLayer(ke_h=0.05, ke_v=0.05, d=3., ks_h = 0.02, ks_v = 0.02)

# S = model.SingleScatRT(surface=soil, canopy=can, models=models, theta=theta, freq=freq)
# S.sigma0()


# f = plt.figure()
# ax = f.add_subplot(111)
# ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0g']['vv']), color='red', linestyle='-', label='vv s0g')
# ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0c']['vv']), color='blue', linestyle='-', label='vv s0c')
# ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0cgt']['vv']), color='green', linestyle='-', label='vv s0cgt')
# ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['stot']['vv']), color='black', linestyle='-', label='vv stot')
# ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0gcg']['vv']), color='black', linestyle='--', label='vv s0gcg')
# # ax.plot(np.rad2deg(self.theta), self.h, color='blue', linestyle='--', label='H')
# ax.set_ylim(-50.,0.)
# ax.grid()
# ax.legend(loc = 3)
# plt.show()


theta = np.deg2rad(np.arange(5.,80.))
# theta = np.deg2rad(25.)
# freq = 5.

# stype = 'turbid_rayleigh'
# # stype='turbid_isotropic'
# models = {'surface': 'Oh92', 'canopy': stype}
# eps = 30.8 - 0.j
# d = 3.
# ke = 1.
# omega=0.1
# soil = Soil(eps=eps, f=1., s=0.05)
# can = OneLayer(ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke)
# S = model.SingleScatRT(surface=soil, canopy=can, models=models, theta=theta, freq=freq)
# S.sigma0()

freq = 1.

stype = 'turbid_rayleigh'
#stype='turbid_isotropic'
models = {'surface': 'Oh92', 'canopy': stype}
eps = 15. -3.j
mv = 0.2
d = 0.22
ke = 1.
omega=0.2
soil = Soil(eps=eps, f=1., s=0.02)
can = OneLayer(ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke)
S = model.SingleScatRT(surface=soil, canopy=can, models=models, theta=theta, freq=freq)
S.sigma0()


f = plt.figure()
ax = f.add_subplot(111)
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0g']['vv']), color='red', linestyle='-', label='vv s0g')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0c']['vv']), color='blue', linestyle='-', label='vv s0c')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0cgt']['vv']), color='green', linestyle='-', label='vv s0cgt')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['stot']['vv']), color='black', linestyle='-', label='vv stot')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0gcg']['vv']), color='black', linestyle='--', label='vv s0gcg')
#ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['stot']['hh']), color='orange', linestyle='-', label='hh stot')
#ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['stot']['hv']), color='purple', linestyle='-', label='hv stot')
# ax.plot(np.rad2deg(self.theta), self.h, color='blue', linestyle='--', label='H')
#ax.set_ylim(-14.,-5.)
ax.set_xlim(0.,70.)
ax.grid()
ax.legend(loc = 3)
plt.show()

# print "yes"


# xxx=Prevot93_surface(mv=0.25, theta=np.radians(23), C2=0.153, D=0.304, C1=-11.2)


# yyy=Prevot93_vegetation(theta=np.radians(23), V=2, A=0.056, B=0.423)

# a = xxx.__dict__['surface'] * yyy.__dict__['tau'] + yyy.__dict__['vegetation']




