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



import pylab
import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator

import matplotlib.ticker


import numpy as np

import pdb

from sense import model
import matplotlib.pyplot as plt
from sense.soil import Soil
from sense.canopy import OneLayer
import pandas as pd
import os
from sense.canopy import OneLayer

from scipy.optimize import minimize


# Import data

path = '/media/tweiss/Daten'
name = 'norm_multi'

df = pd.io.parsers.read_csv(os.path.join(path, name + '.csv'), header=[0, 1], sep=';')
df = df.set_index(pd.to_datetime(df['515']['date']))
df = df.drop(df.filter(like='date'), axis=1)

field_508 = df.filter(like='508')

# field_508 = field_508[[(check == 44 or check == 117) for check in field_508[('508','relativeorbit')]]]


lai_508 = field_508.filter(like='LAI')
field_508 = field_508[~np.isnan(lai_508.values)]

# field_508 = field_508[0:46]
# field_508 = field_508[46:77]

theta_508 = field_508.filter(like='theta')
sm_508 = field_508.filter(like='SM')
height_508 = field_508.filter(like='Height')
vwc_508 = field_508.filter(like='VWC')

lai_508 = field_508.filter(like='LAI')
lai_508 = lai_508[~np.isnan(lai_508.values)]
# lai_508 = vwc_508

vv = field_508.filter(like='sigma_sentinel_vv')
vv = vv[~np.isnan(lai_508.values)]
vv = vv.values.flatten()
vh = field_508.filter(like='sigma_sentinel_vh')
vh = vh[~np.isnan(lai_508.values)]
vh = vh.values.flatten()

theta = np.deg2rad(theta_508.values.flatten())
freq = 5.


# eps = 15. -0.j
s = 0.015
# d = 0.22
d = height_508.values.flatten()/100
coef = 0.9
omega = 0.045

# A_hh = -0.46323766
# B_hh = -0.07569564
# A_vv = 0.01
# B_vv = 0.03623961
# C_hh = -13.19637386
# C_vv = -19.99992029
# D_hh = 14.01814786
# D_vv = 11.06995661

# A_vv = 3.40631410e-03
# B_vv = 2.00000000e+00
# C_vv = 1.99999999e+01
# D_vv = 1.00000000e+01

A_hh = -0.46323766
B_hh = -0.07569564
A_vv = -0.43408517
B_vv = -0.04186564
C_hh = -13.19637386
C_vv = -13.18550537
D_hh = 14.01814786
D_vv = 14.07248098


# 508 vh -0.44968113  -0.09535889 -13.20329416  13.98352766


# A_hh = -0.44968113
# B_hh = -0.09535889
# A_vv = -0.41514708
# B_vv = -0.02384012
# C_hh = -13.20329416
# C_vv = -13.17912135
# D_hh = 13.98352766
# D_vv = 14.10434349

# 301 -0.41514708  -0.02384012 -13.17912135  14.10434349

# 508 -0.44457973  -0.05681031 -13.18990776  14.05050507



# A_vv = 0.020
# B_vv = 0.20

# # models = {'surface': 'WaterCloud', 'canopy': 'water_cloud'}
# # models = {'surface': 'Oh92', 'canopy': 'water_cloud'}
# soil = Soil(mv=sm_508.values.flatten(), C_hh=C_hh, C_vv=C_vv, D_hh=D_hh, D_vv=D_vv, V2=lai_508.values.flatten(), s=s, clay=0.3, sand=0.4, f=freq)
# # soil = Soil(mv=sm_508.values.flatten(), f=freq, s=s, clay=0.11, sand=0.11, bulk=1.65)
# can = OneLayer(V1=lai_508.values.flatten(), V2=lai_508.values.flatten(), A_hh=A_hh, B_hh=B_hh, A_vv=A_vv, B_vv=B_vv)
# S_water = model.WaterCloud(surface=soil, canopy=can, models=models, theta=theta, freq=freq)
# S_water.sigma0()


# def solve_fun(coef):

#     # stype = 'turbid_rayleigh'
#     stype='turbid_isotropic'
#     models = {'surface': 'Oh92', 'canopy': stype}

#     ke = coef * lai_508.values.flatten()
#     ke = coef * vwc_508.values.flatten()


#     # soil = Soil(eps=eps, f=freq, s=s)
#     soil = Soil(mv=sm_508.values.flatten(), f=freq, s=s, clay=0.11, sand=0.11, bulk=1.65)
#     can = OneLayer(ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke)

#     S = model.SingleScatRT(surface=soil, canopy=can, models=models, theta=theta, freq=freq)

#     S.sigma0()

#     return S.__dict__['stot']['hv'][0]

# def fun_opt(VALS):
#     return(np.sum(np.square(solve_fun(VALS[0])-vh)))

# guess = [3.9]

# res = minimize(fun_opt,guess,bounds=[(0.00001,200.)])


def solve_fun(coef):

    # stype = 'turbid_rayleigh'
    stype='turbid_isotropic'
    models = {'surface': 'Oh92', 'canopy': stype}

    ke = coef * lai_508.values.flatten()


    # soil = Soil(eps=eps, f=freq, s=s)
    soil = Soil(mv=sm_508.values.flatten(), f=freq, s=s, clay=0.3, sand=0.4, bulk=1.65)
    can = OneLayer(ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke)

    S = model.SingleScatRT(surface=soil, canopy=can, models=models, theta=theta, freq=freq)

    S.sigma0()

    return S.__dict__['stot']['vv'][0]

def fun_opt(VALS):
    # print(VALS[0])
    # print(solve_fun(VALS[0])-vv)
    return(np.sum(np.square(solve_fun(VALS[0])-vv)))

guess = [0.1]


lai_508_old = lai_508
sm_508_old = sm_508
vv_old = vv
vh_old = vh
theta_old = theta
d_old = d

aaa = []

n = 3

for i in range(len(lai_508_old)-n+1):
    lai_508 = lai_508_old[i:i+n]
    sm_508 = sm_508_old[i:i+n]
    vv = vv_old[i:i+n]
    vh = vh_old[i:i+n]
    theta = theta_old[i:i+n]
    d = d_old[i:i+n]
    res = minimize(fun_opt,guess,bounds=[(0.0001,200.)])
    fun_opt(res.x)
    aaa.append(res.x[0])



bbb = np.asarray(aaa)
print(bbb.mean())
print(bbb.std())
print(bbb.min())
print(bbb.max())

pdb.set_trace()
xxx = solve_fun(res.x[0], res.x[1])
xxx_508 = solve_fun(1.10240852,  0.04536135)
xxx_508_vh = solve_fun(0.30488159,  0.0116789)
xxx_301 = solve_fun(0.5838054 ,  0.07593245)
# res=solve_fun(coef=coef,omega=omega)

pdb.set_trace()
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

stot = S.__dict__['stot']['hh'][0][~np.isnan(S.__dict__['stot']['hh'][0])]
plt.plot(field_508.index[~np.isnan(S.__dict__['stot']['hh'][0])], stot)
# plt.plot(S.__dict__['s0g']['hh'])
# plt.plot(S.__dict__['s0c']['vv'])
plt.plot(field_508.index[~np.isnan(S.__dict__['stot']['hh'][0])],10*np.log10(vh[~np.isnan(S.__dict__['stot']['hh'][0])]))
plt.show()


stot = S.__dict__['stot']['vv'][0][~np.isnan(S.__dict__['stot']['vv'][0])]
plt.plot(field_508.index[~np.isnan(S.__dict__['stot']['vv'][0])], stot)
# plt.plot(S.__dict__['s0g']['hh'])
# plt.plot(S.__dict__['s0c']['vv'])
plt.plot(field_508.index[~np.isnan(S.__dict__['stot']['vv'][0])],10*np.log10(vv[~np.isnan(S.__dict__['stot']['vv'][0])]))
plt.show()










stot = 10*np.log10(S.__dict__['stot']['hh'][0][~np.isnan(S.__dict__['stot']['hh'][0])])
plt.plot(field_508.index[~np.isnan(S.__dict__['stot']['hh'][0])], stot)
# plt.plot(S.__dict__['s0g']['hh'])
# plt.plot(S.__dict__['s0c']['vv'])
plt.plot(field_508.index[~np.isnan(S.__dict__['stot']['hh'][0])],10*np.log10(vh[~np.isnan(S.__dict__['stot']['hh'][0])]))
plt.show()


stot = 10*np.log10(S.__dict__['stot']['vv'][0][~np.isnan(S.__dict__['stot']['vv'][0])])
plt.plot(field_508.index[~np.isnan(S.__dict__['stot']['vv'][0])], stot)
# plt.plot(S.__dict__['s0g']['hh'])
# plt.plot(S.__dict__['s0c']['vv'])
plt.plot(field_508.index[~np.isnan(S.__dict__['stot']['vv'][0])],10*np.log10(vv[~np.isnan(S.__dict__['stot']['vv'][0])]))
plt.show()

stot = 10*np.log10(xxx)
plt.plot(field_508.index, stot)
# plt.plot(S.__dict__['s0g']['hh'])
# plt.plot(S.__dict__['s0c']['vv'])
plt.plot(field_508.index,10*np.log10(vv))
plt.show()


stot = 10*np.log10(xxx)
plt.plot(field_508.index, stot)
# plt.plot(S.__dict__['s0g']['hh'])
# plt.plot(S.__dict__['s0c']['vv'])
plt.plot(field_508.index,10*np.log10(vh))
plt.show()





field_508 = df.filter(like='301')

lai_508 = field_508.filter(like='LAI')
field_508 = field_508[~np.isnan(lai_508.values)]

# field_508 = field_508[0:46]
# field_508 = field_508[46:77]

theta_508 = field_508.filter(like='theta')
sm_508 = field_508.filter(like='SM')
height_508 = field_508.filter(like='Height')
lai_508 = field_508.filter(like='LAI')



vv_301 = field_508.filter(like='sigma_sentinel_vv')
vv_301 = vv_301[~np.isnan(lai_508.values)]
vv_301 = vv_301.values.flatten()
vh_301 = field_508.filter(like='sigma_sentinel_vh')
vh_301 = vh_301[~np.isnan(lai_508.values)]
vh_301 = vh_301.values.flatten()


stot = 10*np.log10(xxx_508)
plt.plot(field_508.index, stot)
plt.plot(field_508.index,10*np.log10(vv))
stot = 10*np.log10(xxx_301)
plt.plot(field_508.index, stot)
plt.plot(field_508.index,10*np.log10(vv_301))

plt.show()




stot = 10*np.log10(xxx)
plt.plot(field_508.index, stot)
plt.plot(field_508.index,10*np.log10(vh))
stot = 10*np.log10(xxx_301)
plt.plot(field_508.index, stot)
plt.plot(field_508.index,10*np.log10(vh_508))












plt.plot(field_508.index[~np.isnan(S.__dict__['stot']['vv'][0])], stot)



















fig, ax = plt.subplots(figsize=(20, 10))

plt.ylabel('Backscatter [dB]', fontsize=15)
plt.tick_params(labelsize=12)

ax.plot(lai_508.index, 10*np.log10(vv), 'ks-', label='Sentinel-1 Data Polarisation VV', linewidth=3)
ax.plot(lai_508.index, 10*np.log10(xxx_301), 'bs-', label='Oh92+SSRT')
ax.plot(lai_508.index, S_water.__dict__['stot']['vv'][0], 'rs-', label='Water Cloud Model')
ax.legend()
ax.legend(loc=2, fontsize=12)

ax2 = ax.twinx()
ax2.tick_params(labelsize=12)
ax2.plot(lai_508.index, lai_508.values.flatten(), color='green', label='LAI')
ax2.legend(loc=1, fontsize=12)
ax2.set_ylabel('LAI [m$^3$/m$^3$]', fontsize=15)

ax.grid(linestyle='-', linewidth=1)
ax.grid(b=True, which='minor', linestyle='--', linewidth=0.5)

days = mdates.DayLocator()
ax.xaxis.set_minor_locator(days)

months = MonthLocator()
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b %Y"))

ax.set_ylim([-20,-6])
ax2.set_ylim([0,7])
# ax3.set_xlim(['2017-03-20', '2017-08-08'])

pdb.set_trace()
plt.savefig('/media/tweiss/Work/EGU_plot3_301')






















fig, ax = plt.subplots(figsize=(20, 10))

plt.ylabel('Backscatter [dB]', fontsize=15)
plt.tick_params(labelsize=12)

ax.plot(lai_508.index, 10*np.log10(vh), 'ks-', label='Sentinel-1 Data Polarisation VH', linewidth=3)
ax.plot(lai_508.index, 10*np.log10(xxx_508_vh), 'bs-', label='Oh92+SSRT')
ax.plot(lai_508.index, S_water.__dict__['stot']['hh'][0], 'rs-', label='Water Cloud Model')
ax.legend()
ax.legend(loc=2, fontsize=12)

ax2 = ax.twinx()
ax2.plot(lai_508.index, lai_508.values.flatten(), color='green', label='LAI')
ax2.legend(loc=1, fontsize=12)
ax2.tick_params(labelsize=12)
ax2.set_ylabel('LAI [m$^3$/m$^3$]', fontsize=15)

ax.grid(linestyle='-', linewidth=1)
ax.grid(b=True, which='minor', linestyle='--', linewidth=0.5)

days = mdates.DayLocator()
ax.xaxis.set_minor_locator(days)

months = MonthLocator()
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b %Y"))

ax.set_ylim([-26,-10])
ax2.set_ylim([0,8])
# ax3.set_xlim(['2017-03-20', '2017-08-08'])

pdb.set_trace()
plt.savefig('/media/tweiss/Work/EGU_plot3_508_vh')
