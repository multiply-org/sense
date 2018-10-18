
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
import scipy.stats
from scipy.optimize import minimize
import pdb


# Helper functions for statistical parameters
#--------------------------------------------
def rmse(predictions, targets):
    """ calculation of RMSE """
    return np.sqrt(((predictions - targets) ** 2).mean())

def linregress(predictions, targets):
    """ Calculate a linear least-squares regression for two sets of measurements """
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(predictions, targets)
    return slope, intercept, r_value, p_value, std_err

def read_mni_data(path, file_name, extention, field, sep=';'):
    """ read MNI campaign data """
    df = pd.io.parsers.read_csv(os.path.join(path, file_name + extension), header=[0, 1], sep=sep)
    df = df.set_index(pd.to_datetime(df[field]['date']))
    df = df.drop(df.filter(like='date'), axis=1)
    return df

def filter_relativorbit(data, field, orbit1, orbit2=None, orbit3=None, orbit4=None):
    """ data filter for relativ orbits """
    output = data[[(check == orbit1 or check == orbit2 or check == orbit3 or check == orbit4) for check in data[(field,'relativeorbit')]]]

### Data preparation ###
#-----------------------------------------------------------------
# storage information
path = '/media/tweiss/Daten'
file_name = 'multi' # theta needs to be changed to for norm multi
extension = '.csv'

field = '508'
pol = 'vv'

# Read MNI data
df = read_mni_data(path, file_name, extension, field)

# filter for field
field_data = df.filter(like=field)

# filter for relativorbit
field_data_orbit = filter_relativorbit(field_data, field, 117)

# get rid of NaN values
parameter_nan = 'LAI'
field_data = field_data[~np.isnan(field_data.filter(like=parameter_nan).values)]

# available auxiliary data
theta_field = np.deg2rad(field_data.filter(like='theta'))
# theta_field[:] = 45
sm_field = field_data.filter(like='SM')
height_field = field_data.filter(like='Height')/100
lai_field = field_data.filter(like='LAI')
vwc_field = field_data.filter(like='VWC')
pol_field = field_data.filter(like='sigma_sentinel_'+pol)
#-----------------------------------------------------------------

### Settings SenSe module ###
#-----------------------------------------------------------------
## Choose models
#---------------
# surface = 'Oh92'
# surface = 'Oh04'
# surface = 'Dubois95'
surface = 'WaterCloud'
# surface = 'I2EM'
# canopy = 'turbid_isotropic'
# canopy = 'turbid_rayleigh'
canopy = 'water_cloud'

models = {'surface': surface, 'canopy': canopy}

## Parameters
#----------------
freq = 5.405

### Surface
#----------
#### Water Cloud
#----------------
C_hh = 0
D_hh = 0
C_hv = -13.19637386
D_hv = 14.01814786
C_vv = -13.18550537
D_vv = 14.07248098
V2 = lai_field.values.flatten()

### Canopy
#----------
#### Water Cloud
#----------------
A_hh = 0
B_hh = 0
A_hv = -0.46323766
B_hv = -0.07569564
A_vv = 0.0029
B_vv = 0.33
V1 = lai_field.values.flatten()
V2 = V2 # initialize in surface model

### Optimization ###
#-----------------------------------------------------------------
# def solve_fun(VALS):

#     for i in range(len(var_opt)):
#         dic[var_opt[i]] = VALS[i]

#     # surface
#     soil = Soil(surface=surface, mv=dic['mv'], C_hh=dic['C_hh'], C_vv=dic['C_vv'], D_hh=dic['D_hh'], D_vv=dic['D_vv'], C_hv=dic['C_hv'], D_hv=dic['D_hv'], V2=dic['V2'])

#     # canopy
#     can = OneLayer(canopy=dic['canopy'], V1=dic['V1'], V2=dic['V2'], A_hh=dic['A_hh'], B_hh=dic['B_hh'], A_vv=dic['A_vv'], B_vv=dic['B_vv'], A_hv=dic['A_hv'], B_hv=dic['B_hv'])

#     S = model.RTModel(surface=soil, canopy=can, models=models, theta=dic['theta'], freq=dic['f'])
#     S.sigma0()

#     return S.__dict__['stot'][pol[::-1]]

# def fun_opt(VALS):
#     return(np.sum(np.square(solve_fun(VALS)-dic['pol_value'])))

# n = 9
# aaa = []
# for i in range(len(pol_field.values.flatten())-n+1):

#     dic = {"mv":sm_field.values.flatten()[i:i+n], "C_hh":C_hh, "C_vv":C_vv, "D_hh":D_hh, "D_vv":D_vv, "C_hv":C_hv, "D_hv":D_hv, "V2":V2[i:i+n], "canopy":canopy, "V1":V1[i:i+n], "A_hh":A_hh, "B_hh":B_hh, "A_vv":A_vv, "B_vv":B_vv, "A_hv":A_hv, "B_hv":B_hv, "lai":lai_field.values.flatten()[i:i+n], "vwc":vwc_field.values.flatten()[i:i+n], "pol_value":pol_field.values.flatten()[i:i+n], "theta":theta_field.values.flatten()[i:i+n], "f":freq}

#     var_opt = ['C_vv', 'D_vv']
#     guess = [-13, 15]
#     bounds = [(-200.,200.),(-200.,200.)]

#     # method = 'L-BFGS-B'
#     res = minimize(fun_opt,guess,bounds=bounds)

#     fun_opt(res.x)
#     aaa.append(res.x)
#-----------------------------------------------------------------

### preparation for optimized run
#-----------------------------------------------------------------
# n = np.int(np.floor(n/2))

# field_data = field_data.drop(field_data.index[-n:])
# field_data = field_data.drop(field_data.index[0:n])
# theta_field = theta_field.drop(theta_field.index[-n:])
# theta_field = theta_field.drop(theta_field.index[0:n])

# sm_field = field_data.filter(like='SM')
# height_field = field_data.filter(like='Height')/100
# lai_field = field_data.filter(like='LAI')
# vwc_field = field_data.filter(like='VWC')

# vv_field = field_data.filter(like='sigma_sentinel_vv')
# vh_field = field_data.filter(like='sigma_sentinel_vh')

# pol_field = field_data.filter(like='sigma_sentinel_'+pol)

# C_vv = np.asarray([el[0] for el in aaa])
# D_vv = np.asarray([el[1] for el in aaa])

#-----------------------------------------------------------------

### Run SenSe module
#-----------------------------------------------------------------

soil = Soil(surface=surface, mv=sm_field.values.flatten(), C_hh=C_hh, C_vv=C_vv, D_hh=D_hh, D_vv=D_vv, C_hv=C_hv, D_hv=D_hv, V2=lai_field.values.flatten())

can = OneLayer(canopy=canopy, V1=lai_field.values.flatten(), V2=lai_field.values.flatten(), A_hh=A_hh, B_hh=B_hh, A_vv=A_vv, B_vv=B_vv, A_hv=A_hv, B_hv=B_hv)

S = model.RTModel(surface=soil, canopy=can, models=models, theta=theta_field.values.flatten(), freq=freq)
S.sigma0()
#-----------------------------------------------------------------

### Plotting
#-----------------------------------------------------------------

# Scatterplot
#------------
plt.plot(10*np.log10(pol_field.values.flatten()),10*np.log10(S.__dict__['stot'][pol[::-1]]), 'ks')

plt.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==44]),10*np.log10(S.__dict__['stot'][pol[::-1]][field_data[(field,'relativeorbit')]==44]), 'ys', label=44)
plt.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==117]),10*np.log10(S.__dict__['stot'][pol[::-1]][field_data[(field,'relativeorbit')]==117]), 'ms', label=117)
plt.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==95]),10*np.log10(S.__dict__['stot'][pol[::-1]][field_data[(field,'relativeorbit')]==95]), 'rs', label=95)
plt.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==168]),10*np.log10(S.__dict__['stot'][pol[::-1]][field_data[(field,'relativeorbit')]==168]), 'gs', label=168)
plt.legend()

x = np.linspace(np.min(10*np.log10(pol_field.values.flatten()))-2, np.max(10*np.log10(pol_field.values.flatten()))+2, 16)
plt.plot(x,x)
plt.savefig('/media/tweiss/Daten/plots/scatterplot_'+field+'_'+pol+'_'+file_name+'_'+S.models['surface']+'_'+S.models['canopy'])
plt.close()


# Plot
#------
date = field_data.index

fig, ax = plt.subplots(figsize=(20, 10))
# plt.title('Winter Wheat')
plt.ylabel('Backscatter [dB]', fontsize=15)
plt.tick_params(labelsize=12)

ax.plot(10*np.log10(pol_field), 'ks-', label='Sentinel-1 Pol: ' + pol, linewidth=3)

ax.plot(date, 10*np.log10(S.__dict__['s0g'][pol[::-1]]), 'rs-', label=pol+' s0g')
ax.plot(date, 10*np.log10(S.__dict__['s0c'][pol[::-1]]), 'cs-', label=pol+' s0c')

ax.plot(date, 10*np.log10(S.__dict__['stot'][pol[::-1]]), 'C1s-', label=S.models['surface']+ ' + ' +  S.models['canopy'] + ' Pol: ' + pol)
ax.legend()
ax.legend(loc=2, fontsize=12)


ax6 = ax.twinx()
ax6.tick_params(labelsize=12)
lns1 = ax6.plot(vwc_field, 'g', label='VWC')
ax6.set_ylabel('VWC [kg/m2]', fontsize=15, color='green')
ax6.yaxis.label.set_color('green')
ax6.spines['right'].set_position(('outward', 60))

ax.grid(linestyle='-', linewidth=1)
ax.grid(b=True, which='minor', linestyle='--', linewidth=0.5)

days = mdates.DayLocator()
ax.xaxis.set_minor_locator(days)

months = MonthLocator()
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b %Y"))


ax.set_ylim([-30,-5])

ax.set_xlim(['2017-03-20', '2017-07-15'])


slope, intercept, r_value, p_value, std_err = scipy.stats.linregress((pol_field.values.flatten()), (S.__dict__['stot'][pol[::-1]]))
slope1, intercept1, r_value1, p_value1, std_err1 = scipy.stats.linregress(10*np.log10(pol_field.values.flatten()), 10*np.log10(S.__dict__['stot'][pol[::-1]]))
rmse = rmse(10*np.log10(pol_field.values.flatten()), 10*np.log10(S.__dict__['stot'][pol[::-1]]))

plt.title('Winter Wheat, R2 = ' + str(r_value) + ' RMSE = ' + str(rmse))
plt.savefig('/media/tweiss/Daten/plots/plot_'+field+'_'+pol+'_'+file_name+'_'+S.models['surface']+'_'+S.models['canopy'])
plt.close()



pdb.set_trace()

