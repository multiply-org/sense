"""
compare results of implemented models
against Online Code 10.1 (http://mrs.eecs.umich.edu/codes/Module10_1/Module10_1.html)

real part: 5.46
imaag part: 0.37
frequenz: 1.3
rms: 0.011
correlation length: 0.373
correlation function: exponential or gaussian!!!

I2EM exp15: hh and vv looks pretty good, HV a little bit of
I2EM gauss: hh and vv looks pretty good, HV better as exp15 but still not perfect
"""

import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__)) + os.sep + '..')

import numpy as np

from sense.surface import Dubois95, Oh92, I2EM
from sense.util import f2lam

import matplotlib.pyplot as plt
import pdb
plt.close('all')

theta = np.deg2rad(np.arange(20.,70.))
eps = 5.46 -0.37j
eps = 5.1 -1j

s_vv_i = []
s_hh_i = []
s_hv_i = []

f  = 4.5  # GHz
# lam = f2lam(f)  # m
# k = 2.*np.pi/lam

# ks = 0.4
# s = ks/k
# l=6.*s

s = 0.01
l = 0.041

for THETA in theta:

    theta_one_value = THETA

    # IEM
    acf_type='gauss'
    acf_type='exp15'
    I = I2EM(f, eps, s, l, theta_one_value, acf_type=acf_type, xpol=True)
    s_vv_i.append(I.vv)
    s_hh_i.append(I.hh)
    s_hv_i.append(I.hv)

s_vv_i = 10.*np.log10(np.array(s_vv_i))
s_hh_i = 10.*np.log10(np.array(s_hh_i))
s_hv_i = 10.*np.log10(np.array(s_hv_i))

f = plt.figure()

theta = np.rad2deg(theta)

plt.plot(theta, s_vv_i, color='red', label='VV I2EM')
plt.plot(theta, s_hh_i, color='blue', label='HH I2EM')
plt.plot(theta, s_hv_i, color='green', label='HV I2EM')
plt.grid()
plt.legend()
plt.xlim(0.,70.)
plt.ylim(-50.,30.)
# plt.ylim(-200.,25.)
plt.xlabel('Incidence Angle [°]')
plt.ylabel('backscattering coefficient VV [dB]')





pdb.set_trace()



plt.show()
