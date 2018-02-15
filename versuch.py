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
from sense.canopy import OneLayer, WaterCloudLayer
import pdb


V1=np.asarray([ 0.4211875 ,  0.7969375 ,  1.037     ,  1.41933333,  1.45933333,
        1.51933333,  1.65933333,  1.68383333,  1.63223333,  1.59926667,
        1.54756061,  1.68465152,  1.89028788,  2.37010606,  2.50719697,
        2.71283333,  2.84421212,  3.05698864,  3.08180682,  3.11903409,
        3.20589773,  3.23071591,  3.26794318,  3.29172727,  3.36486979,
        3.39424479,  3.43830729,  3.57288194,  3.75471528,  4.02746528,
        4.20172222,  4.65491667,  4.79371667,  5.00191667,  5.48771667,
        5.62651667,  5.83471667,  6.01768056,  6.01625   ,  5.99645   ,
        5.95025   ,  5.93705   ,  5.91725   ,  5.80010985,  5.72820076,
        5.62033712,  5.3686553 ,  5.29674621,  5.18888258,  5.17076786,
        5.204625  ,  5.25541071,  5.175675  ,  5.045475  ,  4.850175  ,
        4.394475  ,  4.264275  ,  4.068975  ,  4.031     ,  4.031     ,
        4.031     ,  4.031     ,  4.031     ,  4.031     ])
# V1[V1 > 5.] = 0.8
# V1[V1 < 1.] = 0.1
# V1[V1 > 5.] = 0.8
# V1[V1 > 4.] = 0.7
# V1[V1 > 3.] = 0.5
# V1[V1 > 2.] = 0.35
# V1[V1 > 1.] = 0.2

V2=V1

mv=np.asarray([ 0.20838876,  0.2015811 ,  0.20062526,  0.20054675,  0.20092561,
        0.19796436,  0.20694698,  0.20858946,  0.20522828,  0.20463417,
        0.20530814,  0.20213852,  0.20231781,  0.21046817,  0.2391101 ,
        0.28449486,  0.27608734,  0.23792895,  0.24135696,  0.22690808,
        0.3256423 ,  0.31161885,  0.27498145,  0.27151334,  0.27793293,
        0.27486545,  0.26108529,  0.29115406,  0.27638248,  0.25717182,
        0.26358998,  0.29552464,  0.27529122,  0.25019986,  0.26082066,
        0.26009207,  0.25267054,  0.25047114,  0.24810711,  0.23815195,
        0.24828169,  0.29548717,  0.26351374,  0.26228643,  0.25721306,
        0.24093643,  0.24248841,  0.23850948,  0.22759338,  0.2277716 ,
        0.22636968,  0.21697637,  0.22371928,  0.24591281,  0.28826846,
        0.24526919,  0.23444086,  0.21395943,  0.21274064,  0.22877096,
        0.21838316,  0.24939461,  0.23055565,  0.21548489])

# theta=np.deg2rad(np.asarray([ 44.69104004,  42.44337463,  33.25333405,  35.97608948,
#         44.69067764,  42.44624329,  35.99119186,  44.69095612,
#         42.44985199,  33.26148605,  35.97618866,  44.68642426,
#         42.44688797,  35.98065567,  44.69062805,  42.45200729,
#         33.26424789,  35.97745895,  44.68911362,  42.44366074,
#         35.9924469 ,  44.68640518,  42.44782257,  33.26057053,
#         35.98235321,  44.68797684,  42.45107269,  35.98369598,
#         44.6839447 ,  42.44587708,  33.25891876,  35.98817444,
#         44.68839645,  42.45076752,  35.98382568,  44.68654251,
#         42.44100189,  35.99754715,  44.69528198,  42.44817734,
#         35.98278809,  44.69725037,  42.44649887,  35.97967148,
#         44.68753433,  42.44492722,  35.98325348,  44.6935997 ,
#         42.44932175,  35.98846436,  44.68527603,  42.44314194,
#         35.9786377 ,  44.69162369,  42.44630051,  35.98388672,
#         44.69190979,  42.44133377,  35.98124313,  44.6942482 ,
#         42.44450378,  35.98654556,  44.69247055,  42.44791031]))
theta=np.deg2rad(35)

# theta = np.deg2rad(np.arange(1.,80.))
# theta = np.deg2rad(23.)
freq = 5.

# stype = 'turbid_rayleigh'
stype='turbid_isotropic'
# stype = 'water_cloud'
models = {'surface': 'Oh92', 'canopy': stype}
# models = {'surface': 'WaterCloud', 'canopy': stype}


# eps = 15. -0.j
# mv = np.asarray([0.2, 0.3, 0.4])
# mv = 0.2
# eps = [15. -0.j, 15. -.0j]
s=0.02
d = 0.22
ke = 1.
omega=0.2

# V1 = 7.
# V2 = 7.
A_hh = 0.02
B_hh = 0.21
A_vv = 0.01
B_vv = 0.03623961
C_hh = -14.8465
C_vv = -19.99992029
D_hh = 15.907
D_vv = 11.06995661

A_vv = 3.40631410e-03
B_vv = 2.00000000e+00
C_vv = 1.99999999e+01
D_vv = 1.00000000e+01

# soil = Soil(eps=eps, f=freq, s=s)
soil = Soil(mv=mv, f=freq, s=s, clay=0.3, sand=0.4)
# soil = Soil(mv=mv, C_hh=C_hh, C_vv=C_vv, D_hh=D_hh, D_vv=D_vv, V2=V2, s=s, clay=0.3, sand=0.4, f=freq)
can = OneLayer(ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke)
# can = WaterCloudLayer(V1=V1, V2=V2, A_hh=A_hh, B_hh=B_hh, A_vv=A_vv, B_vv=B_vv)
S = model.SingleScatRT(surface=soil, canopy=can, models=models, theta=theta, freq=freq)
# S = model.WaterCloud(surface=soil, canopy=can, models=models, theta=theta, freq=freq)
S.sigma0()
pdb.set_trace()


f = plt.figure()
ax = f.add_subplot(111)
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0g']['vv']), color='red', linestyle='-', label='vv s0g')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0c']['vv']), color='blue', linestyle='-', label='vv s0c')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0cgt']['vv']), color='green', linestyle='-', label='vv s0cgt')
# ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['stot']['vv']), color='black', linestyle='-', label='vv stot')
ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['s0gcg']['vv']), color='green', linestyle='--', label='vv s0gcg')
# ax.plot(np.rad2deg(S.__dict__['theta']), 10*np.log10(S.__dict__['stot']['hh']), color='orange', linestyle='--', label='hh stot')
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


