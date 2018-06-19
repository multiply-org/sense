"""
comparison with figure 10-10
in Ulaby 2014

Comments Alex
------------
- cross pol for Gaussian seems rather low

number of ittertaions??h


in general something with the Gaussian seems not to work yet!

quite o.k. for hv


Comments Thomas
-------------------
hh and vv polarisation for Exponential and Gaussian are the same for python implementation and Ulaby online Code 10.1. Figure 10-10 has slightly different values especially for Gaussian.

Online code for hv polarisation doesn't work except by changing IEMX_model lines:

svh = dblquad(@(r,phi)xpol_integralfunc(r, phi, sp,xx, ks2, cs,s, kl2, L, er, rss, rvh, n_spec), 0.1, 1, 0, pi);

to

svh = dblquad(@(r,phi)xpol_integralfunc(r, phi, sp,xx, ks2, cs,s, kl2, L, er, rss, rvh, n_spec), 0.1, 1.0, 0, 1, 0.5);

(Problem absolute tolerance)

and

VH = 4 * acc .* Fvh .* vhmnsum .*r

to

VH = 4 * acc * Fvh * vhmnsum *r;

(Problem Matrix multiplication)

"""

import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__)) + os.sep + '..')

import numpy as np

from sense.surface import I2EM

import matplotlib.pyplot as plt

import time

import pdb

def db(x):
    return 10.*np.log10(x)

start = time.time()
plt.close('all')


theta_deg = np.linspace(0.,70., 71)
theta = np.deg2rad(theta_deg)

f = plt.figure()
ax = f.add_subplot(111)


eps = 11.3-1.5j
f = 3.

s = 1./100.
l = 10./100.

# eps = 11.3
# theta=np.deg2rad(30)
# xpol = True
# auto=False
# I1 = I2EM(f, eps, s, l, theta, acf_type='gauss', xpol=xpol, auto=auto)
# pdb.set_trace()
# eps = 3-0j
# f = 2
# t = np.deg2rad(60)
# xpol = True
# auto=False

# I1 = I2EM(f, eps, s, l, t, acf_type='exp15', xpol=xpol, auto=auto)
# 10*np.log10(I1._i2em_bistatic())
# pdb.set_trace()

# hh1=[]
# hh2=[]
# vv1=[]
# vv2=[]
# hv1=[]
# hv2=[]
# xpol = True
# auto=False
# for t in theta:
#     I1 = I2EM(f, eps, s, l, t, acf_type='gauss', xpol=xpol, auto=auto)
#     I2 = I2EM(f, eps, s, l, t, acf_type='exp15', xpol=xpol, auto=auto)
#     hh1.append((I1._i2em_bistatic()[0]))
#     vv1.append((I1._i2em_bistatic()[1]))
#     hh2.append((I2._i2em_bistatic()[0]))
#     vv2.append((I2._i2em_bistatic()[1]))



# plt.plot(theta_deg, 10*np.log10(hh1))
# plt.plot(theta_deg, 10*np.log10(hh2))
# plt.plot(theta_deg, 10*np.log10(vv1))
# plt.plot(theta_deg, 10*np.log10(vv2))
# pdb.set_trace()

# theta2 = [10,20,30,40,50,60,70]
# hhh = [-6.9240, -12.733, -16.440, -19.085, -21.353, -23.941, -27.973]
# plt.plot(theta2, hhh)

# hm = -6.9240
# hmno = -7.0843

# hm = -12.733
# hmno = -13.463

# hm = -16.440
# hmno = -18.069

# hm = -19.085
# hmno = -21.770

# hm = -21.353
# hmno = -25.011

# hm = -23.941
# hmno = -28.039

# hm = -27.973
# hmno = -30.923

hh1=[]
hh2=[]
vv1=[]
vv2=[]
hv1=[]
hv2=[]
xpol = True
auto=False
for t in theta:
    print(t)
    I1 = I2EM(f, eps, s, l, t, acf_type='gauss', xpol=xpol, auto=auto)
    I2 = I2EM(f, eps, s, l, t, acf_type='exp15', xpol=xpol, auto=auto)
    print(I1.ks, I1.kl)
    hh1.append(I1.hh)
    hh2.append(I2.hh)
    vv1.append(I1.vv)
    vv2.append(I2.vv)
    if xpol:
        hv1.append(I1.hv)
        hv2.append(I2.hv)

hh1 = np.array(hh1)
hh2 = np.array(hh2)
vv1 = np.array(vv1)
vv2 = np.array(vv2)
hv1 = np.array(hv1)
hv2 = np.array(hv2)

# ax.plot(theta_deg, db(hh2), color='red', label='hh')
ax.plot(theta_deg, db(hh1), color='blue', label='hh')

# ax.plot(theta_deg, db(vv2), color='red', label='vv', linestyle='--')
ax.plot(theta_deg, db(vv1), color='blue', label='vv', linestyle='--')

# ax.plot(theta_deg, db(hv2), color='red', label='hv', linestyle='-.')
ax.plot(theta_deg, db(hv1), color='blue', label='hv', linestyle='-.')
pdb.set_trace()
ax.grid()
ax.set_xlim(0.,70.)
ax.set_ylim(-100.,15.)
ax.set_yticks(np.arange(-90, 20, 10))
print('Elapsed time [s]: ', time.time()-start)
plt.show()
pdb.set_trace()
