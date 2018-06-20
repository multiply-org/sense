
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator
# import matplotlib.ticker
import numpy as np
from sense.canopy import OneLayer
from sense.soil import Soil
from sense import model

from scipy.optimize import minimize


import pdb


# Import auxiliary data
#-----------------------
path = '/media/tweiss/Daten'
file_name = 'multi'
field = '508'
pol = 'vv'
pol2 = 'vv'
"""theta needs to be changed to for norm multi'!!!!!!!!!!!!!!!!"""

# Read auxiliary data
#----------------------
df = pd.io.parsers.read_csv(os.path.join(path, file_name + '.csv'), header=[0, 1], sep=';')
df = df.set_index(pd.to_datetime(df[field]['date']))
df = df.drop(df.filter(like='date'), axis=1)

# filter for field
#-------------------
field_data = df.filter(like=field)

# get rid of NaN values
#------------------------
parameter_nan = 'LAI'
field_data = field_data[~np.isnan(field_data.filter(like=parameter_nan).values)]

# field_data = field_data[field_data[(field,'relativeorbit')]==44]
# field_data = field_data[[(check == 44 or check == 117) for check in field_data[(field,'relativeorbit')]]]

# n = 1
# field_data = field_data.drop(field_data.index[-n:])
# field_data = field_data.drop(field_data.index[0:n])


# available auxiliary data
#--------------------------
theta_field = field_data.filter(like='theta')
# theta_field[field,'theta']=35.
theta_field = np.deg2rad(theta_field)

sm_field = field_data.filter(like='SM')
height_field = field_data.filter(like='Height')/100
lai_field = field_data.filter(like='LAI')
vwc_field = field_data.filter(like='VWC')

# vv_field = field_data.filter(like='sigma_sentinel_vv')
# vv = vv_field.values.flatten()
# vh_field = field_data.filter(like='sigma_sentinel_vh')
# vh = vh_field.values.flatten()

pol_field = field_data.filter(like='sigma_sentinel_'+pol)


# Settings SenSe module
#-----------------------
## Parameters
#----------------
freq = 5.

### Surface
#----------------
#### Water Cloud
#----------------
C_hh = -13.19637386
D_hh = 14.01814786
C_vv = -13.18550537
D_vv = 14.07248098

#### Oh92
#--------
clay = 0.3
sand = 0.4
bulk = 1.65
s = 0.015

#### Dubois
#-----------

### Canopy
#---------
#### Water Cloud
#----------------
A_hh = -0.46323766
B_hh = -0.07569564
A_vv = -0.43408517
B_vv = -0.04186564

#### SSRT canopy
#---------------
coef = 1.5
omega = 0.005
# coef = 1.10240852
# omega = 0.04536135
# coef = np.arange(len(theta_field), dtype=float)
# coef[0:35] = 1.10240852
# coef[35:len(coef)] = 0.3

## Choose models
#---------------
surface = 'Oh92'
surface = 'Oh04'
surface = 'Dubois95'
surface = 'WaterCloud'
canopy = 'turbid_isotropic'
# canopy = 'turbid_rayleigh'
# canopy = 'water_cloud'

models = {'surface': surface, 'canopy': canopy}




# Optimisation
#--------------

# def solve_fun(coef,omega):

#     # stype = 'turbid_rayleigh'
#     stype='turbid_isotropic'
#     models = {'surface': 'Oh92', 'canopy': stype}

#     ke = coef * np.sqrt(lai)
#     ke = coef * np.sqrt(vwc)
#     omega = 0.045

#     # soil = Soil(eps=eps, f=freq, s=s)
#     soil = Soil(mv=sm, f=freq, s=s, clay=0.3, sand=0.4, bulk=1.65)
#     can = OneLayer(ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke)

#     S = model.SingleScatRT(surface=soil, canopy=can, models=models, theta=theta, freq=freq)

#     S.sigma0()

#     return S.__dict__['stot'][pol2][0]

# def fun_opt(VALS):
#     # pdb.set_trace()
#     return(np.sum(np.square(solve_fun(VALS[0],VALS[1])-pol_value)))

# guess = [0.1, 0.045]


def solve_fun_SSRT(coef):

    # ke = coef * np.sqrt(lai)
    ke = coef * np.sqrt(vwc)

    # initialize surface
    #--------------------
    soil = Soil(mv=sm, C_hh=C_hh, C_vv=C_vv, D_hh=D_hh, D_vv=D_vv, V2=lai_field.values.flatten(), s=s, clay=clay, sand=sand, f=freq, bulk=bulk)

    # initialize canopy
    #-------------------
    can = OneLayer(canopy=canopy, ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke, V1=lai_field.values.flatten(), V2=lai_field.values.flatten(), A_hh=A_hh, B_hh=B_hh, A_vv=A_vv, B_vv=B_vv)

    # run SenSe module
    #------------------
    S = model.RTModel(surface=soil, canopy=can, models=models, theta=theta, freq=freq)
    S.sigma0()

    return S.__dict__['stot'][pol2][0]

def fun_opt(VALS):
    return(np.sum(np.square(solve_fun_SSRT(VALS[0])-pol_value)))

guess = [0.45]


# def solve_fun_surfacewatercloud(coef, C_vv, D_vv):

#     # ke = coef * np.sqrt(lai)
#     ke = coef * np.sqrt(vwc)

#     # initialize surface
#     #--------------------
#     soil = Soil(mv=sm, C_hh=C_hh, C_vv=C_vv, D_hh=D_hh, D_vv=D_vv, V2=lai_field.values.flatten(), s=s, clay=clay, sand=sand, f=freq, bulk=bulk)

#     # initialize canopy
#     #-------------------
#     can = OneLayer(canopy=canopy, ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke, V1=lai_field.values.flatten(), V2=lai_field.values.flatten(), A_hh=A_hh, B_hh=B_hh, A_vv=A_vv, B_vv=B_vv)

#     # run SenSe module
#     #------------------
#     S = model.RTModel(surface=soil, canopy=can, models=models, theta=theta, freq=freq)
#     S.sigma0()

#     return S.__dict__['stot'][pol2][0]

# def fun_opt(VALS):
#     return(np.sum(np.square(solve_fun_surfacewatercloud(VALS[0],VALS[1],VALS[2])-pol_value)))

# guess = [0.45, -13., 14.]













lai_508_old = lai_field.values.flatten()
vwc_508_old = vwc_field.values.flatten()
sm_508_old = sm_field.values.flatten()
pol_old = pol_field.values.flatten()
# vv_old = vv_field.values.flatten()
# vh_old = vh_field.values.flatten()
theta_old = theta_field.values.flatten()
d_old = height_field.values.flatten()

aaa = []
bbb = []
ccc = []

n = 3
for i in range(len(lai_508_old)-n+1):
    lai = lai_508_old[i:i+n]
    vwc = vwc_508_old[i:i+n]
    sm = sm_508_old[i:i+n]
    pol_value = pol_old[i:i+n]
    # vv = vv_old[i:i+n]
    # vh = vh_old[i:i+n]
    theta = theta_old[i:i+n]
    d = d_old[i:i+n]
    # res = minimize(fun_opt,guess,bounds=[(0.0001,200.),(0.0001,200.)], method='L-BFGS-B')
    res = minimize(fun_opt,guess,bounds=[(0.001,200.)])
    # res = minimize(fun_opt,guess,bounds=[(0.001,200.),(-100.,200.),(-100.,200.)])
    fun_opt(res.x)
    aaa.append(res.x[0])
    # bbb.append(res.x[1])
    # ccc.append(res.x[2])

n = 1
field_data = field_data.drop(field_data.index[-n:])
field_data = field_data.drop(field_data.index[0:n])

# available auxiliary data
#--------------------------
theta_field = field_data.filter(like='theta')
# theta_field[field,'theta']=10.
theta_field = np.deg2rad(theta_field)

sm_field = field_data.filter(like='SM')
height_field = field_data.filter(like='Height')/100
lai_field = field_data.filter(like='LAI')
vwc_field = field_data.filter(like='VWC')

vv_field = field_data.filter(like='sigma_sentinel_vv')
# vv = vv_field.values.flatten()
vh_field = field_data.filter(like='sigma_sentinel_vh')
# vh = vh_field.values.flatten()

pol_field = field_data.filter(like='sigma_sentinel_'+pol)

coef = aaa
# omega = bbb
# omega= 0.045
# coef = [2.8211794522771836, 3.5425090255370901, 3.2928811420649549, 3.7054736773563386, 2.8961183027242945, 2.6151602913331349, 2.4189863720217848, 2.4167831184418151, 2.9416289237554785, 3.173742070422747, 3.2383792609197575, 2.7931417278929307, 3.1121882936071223, 3.1256520435841559, 3.4783185812129198, 3.221982368854273, 3.3728926537091137, 3.3781506298071911, 3.2569473671053992, 2.4663320815477578, 2.1141754601248004, 2.1209806250929346, 2.0649019533248212, 2.9544567022270387, 3.2224869016201918, 3.2369640163741287, 3.1260764088434825, 2.2886755717903537, 2.444457485088495, 2.4480123008899097, 2.9424868137311013, 2.6918809874287475, 2.7294523688168129, 2.797705352991366, 2.0186456268959168, 2.8744575821617224, 3.044882444109152, 2.9802180637002609, 2.6019756936977609, 2.5503946347668287, 2.3604725553358343, 2.3317510351411883, 1.9275854991156061, 2.2753308530216287, 2.4092428831073858, 2.3698519298054159, 1.9286656987151718, 1.6963728219945839, 2.0918823569585152, 1.9957476741875333, 1.7201838997262218, 1.4695235195381868, 1.5361435193274382, 1.6667889407678251, 0.82853524517534494, 0.79365611796719404, 0.83635432955114541, 0.74755549562741264, 0.52765406146373417, 0.48239419873329636, 0.54280068255034319, 0.5049404990262979, 0.38576228044679073, 0.45119396015499219, 0.5469751828940197, 0.57279234317771832, 0.38171780066089117, 0.34532639205524862, 0.31494432189245669, 0.29321333743882261, 0.17349153645934853, 0.23334377236206902, 0.30119514955623677, 0.35404564503560393]


# ke = coef * np.sqrt(lai_field.values.flatten())
ke = coef * np.sqrt(vwc_field.values.flatten())



# C_vv = bbb
# D_vv = ccc





# initialize surface
#--------------------
soil = Soil(mv=sm_field.values.flatten(), C_hh=C_hh, C_vv=C_vv, D_hh=D_hh, D_vv=D_vv, V2=lai_field.values.flatten(), s=s, clay=clay, sand=sand, f=freq, bulk=bulk)

# initialize canopy
#-------------------
can = OneLayer(canopy=canopy, ke_h=ke, ke_v=ke, d=height_field.values.flatten(), ks_h = omega*ke, ks_v = omega*ke, V1=lai_field.values.flatten(), V2=lai_field.values.flatten(), A_hh=A_hh, B_hh=B_hh, A_vv=A_vv, B_vv=B_vv)

# run SenSe module
#------------------
S = model.RTModel(surface=soil, canopy=can, models=models, theta=theta_field.values.flatten(), freq=freq)
S.sigma0()

pdb.set_trace()

# Scatterplot
#------------

plt.plot(10*np.log10(pol_field.values.flatten()),10*np.log10(S.__dict__['stot'][pol2][0]), 'ks')

plt.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==44]),10*np.log10(S.__dict__['stot'][pol2][0][field_data[(field,'relativeorbit')]==44]), 'ys', label=44)
plt.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==117]),10*np.log10(S.__dict__['stot'][pol2][0][field_data[(field,'relativeorbit')]==117]), 'ms', label=117)
plt.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==95]),10*np.log10(S.__dict__['stot'][pol2][0][field_data[(field,'relativeorbit')]==95]), 'rs', label=95)
plt.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==168]),10*np.log10(S.__dict__['stot'][pol2][0][field_data[(field,'relativeorbit')]==168]), 'gs', label=168)
plt.legend()

x = np.linspace(np.min(10*np.log10(pol_field.values.flatten()))-2, np.max(10*np.log10(pol_field.values.flatten()))+2, 16)
plt.plot(x,x)
plt.savefig('/media/tweiss/Daten/plots/scatterplot_'+field+'_'+pol+'_'+file_name+'_'+S.models['surface']+'_'+S.models['canopy'])
plt.close()




# Plot
#------
date = field_data.index

fig, ax = plt.subplots(figsize=(20, 10))
plt.title('Winter Wheat')
plt.ylabel('Backscatter [dB]', fontsize=15)
plt.tick_params(labelsize=12)

ax.plot(10*np.log10(pol_field), 'ks-', label='Sentinel-1 Pol: ' + pol, linewidth=3)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==44]), 'ys-', label=44)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==117]), 'ms-', label=117)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==95]), 'rs-', label=95)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==168]), 'gs-', label=168)

# ax.plot(date, 10*np.log10(S.__dict__['s0g'][pol2]), 'rs-', label=pol+' s0g')
# ax.plot(date, 10*np.log10(S.__dict__['s0c'][pol2]), 'cs-', label=pol+' s0c')
# ax.plot(date, 10*np.log10(S.__dict__['s0cgt'][pol2]), 'ms-', label=pol+' s0cgt')
# ax.plot(date, 10*np.log10(S.__dict__['s0gcg'][pol2]), 'ys-', label=pol+' s0gcg')
ax.plot(date, 10*np.log10(S.__dict__['stot'][pol2][0]), 'C1s-', label=S.models['surface']+ ' + ' +  S.models['canopy'] + ' Pol: ' + pol)
ax.legend()
ax.legend(loc=2, fontsize=12)

# ax2 = ax.twinx()
# ax2.plot(lai_field, color='green', label='LAI')
# ax2.plot(vwc_field, color='blue', label='VWC')
# ax2.legend(loc=1, fontsize=12)
# ax2.tick_params(labelsize=12)
# ax2.set_ylabel('LAI [m$^3$/m$^3$]', fontsize=15, color='green')
# ax3 = ax2.twinx()
# ax3.set_ylabel('VWC [kg/m$^2$]', fontsize=15, color='blue')

# ax3.plot(date,coef)
# ax3.plot(date,ke)
# ax4 = ax3.twinx()
# ax4.plot(sm_field)
# ax4.set_ylim([-0.8,0.4])


ax.grid(linestyle='-', linewidth=1)
ax.grid(b=True, which='minor', linestyle='--', linewidth=0.5)

days = mdates.DayLocator()
ax.xaxis.set_minor_locator(days)

months = MonthLocator()
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b %Y"))

ax.set_ylim([np.min(10*np.log10(pol_field.values.flatten()))-2, np.max(10*np.log10(pol_field.values.flatten()))+2])
# ax2.get_yaxis().set_ticks([])
# ax2.set_ylim([0,8])
# ax3.set_ylim([0,8])
# ax3.set_xlim(['2017-03-20', '2017-08-08'])


plt.savefig('/media/tweiss/Daten/plots/plot_'+field+'_'+pol+'_'+file_name+'_'+S.models['surface']+'_'+S.models['canopy'])
plt.close()

