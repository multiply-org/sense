
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
from scipy.optimize import minimize, root
import pdb
import datetime
import random
import gp_emulator


# Helper functions for statistical parameters
#--------------------------------------------
def rmse_prediction(predictions, targets):
    """ calculation of RMSE """
    return np.sqrt(np.nanmean((predictions - targets) ** 2))

def linregress(predictions, targets):
    """ Calculate a linear least-squares regression for two sets of measurements """
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(predictions, targets)
    return slope, intercept, r_value, p_value, std_err

def read_mni_data(path, file_name, extention, field, sep=','):
    """ read MNI campaign data """
    df = pd.io.parsers.read_csv(os.path.join(path, file_name + extension), header=[0, 1], sep=sep)
    pd.to_datetime(df[field]['date']).dt.date
    df = df.set_index(pd.to_datetime(df[field]['date']).dt.date)

    dd = pd.io.parsers.read_csv(os.path.join('/media/tweiss/Daten/buffer/params_S2_rfieldbuffer50.csv'), header=[0, 1], sep=',')
    ddd = dd.set_index(pd.to_datetime(dd['508_high']['date']).dt.date)
    merge=pd.merge(df,ddd, how='inner', left_index=True, right_index=True)
    df = merge
    #df = df.set_index(pd.to_datetime(df[field]['date']))
    df = df.drop(df.filter(like='date'), axis=1)
    return df

def read_agrometeo(path, file_name, extentio, sep=';', decimal=','):
    """ read agro-meteorological station (hourly data) """
    df = pd.read_csv(os.path.join(path, file_name + extension), sep=sep, decimal=decimal)
    df['SUM_NN050'] = df['SUM_NN050'].str.replace(',','.')
    df['SUM_NN050'] = df['SUM_NN050'].str.replace('-','0').astype(float)

    df['date'] = df['Tag'] + ' ' + df['Stunde']

    df = df.set_index(pd.to_datetime(df['date'], format='%d.%m.%Y %H:%S'))
    return df

def filter_relativorbit(data, field, orbit1, orbit2=None, orbit3=None, orbit4=None):
    """ data filter for relativ orbits """
    output = data[[(check == orbit1 or check == orbit2 or check == orbit3 or check == orbit4) for check in data[(field,'relativeorbit')]]]
    return output

def smooth(x,window_len=1,window='hanning'):
        if x.ndim != 1:
                raise ValueError #, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError #, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError #, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

def read_data(path, file_name, extension, field, path_agro, file_name_agro, extension_agro):
    # Read MNI data
    df = read_mni_data(path, file_name, extension, field)

    # Read agro-meteorological station
    #df_agro = read_agrometeo(path_agro, file_name_agro, extension_agro)
    df_agro=0

    # filter for field
    field_data = df.filter(like=field)

    # filter for relativorbit
    # field_data_orbit = filter_relativorbit(field_data, field+'_x', 168)
    # field_data = field_data_orbit
    field_data_orbit=0

    # get rid of NaN values
    parameter_nan = 'LAI'
    field_data = field_data[~np.isnan(field_data.filter(like=parameter_nan).values)]

    # available auxiliary data
    theta_field = np.deg2rad(field_data.filter(like='theta'))
    # theta_field[:] = 45
    sm_field = field_data.filter(like='SM')
    height_field = field_data.filter(like='Height')/100
    height_field2 = field_data.filter(like='Cab')/100
    lai_field = field_data.filter(like='LAI')
    lai_field2 = field_data.filter(like='S2')
    vwc_field = field_data.filter(like='VWC')
    pol_field = field_data.filter(like='sigma_sentinel_'+pol)

    return df, df_agro, field_data, field_data_orbit, theta_field, sm_field, height_field, lai_field, vwc_field, pol_field, height_field2, lai_field2


### Data preparation ###
#-----------------------------------------------------------------
# storage information
path = '/media/tweiss/Daten/buffer'
file_name = 'multi_fieldbuffer50' # theta needs to be changed to for norm multi
extension = '.csv'

path_agro = '/media/nas_data/2017_MNI_campaign/field_data/meteodata/agrarmeteorological_station'
file_name_agro = 'Eichenried_01012017_31122017_hourly'
extension_agro = '.csv'

field = '508_high'
pol = 'vh'

df, df_agro, field_data, field_data_orbit, theta_field, sm_field, height_field, lai_field, vwc_field, pol_field, height_field2, lai_field2 = read_data(path, file_name, extension, field, path_agro, file_name_agro, extension_agro)


date = pd.to_datetime(field_data.index)

theta = field_data.filter(like='theta').values.flatten()

theta = np.deg2rad(theta)

LAI = field_data.filter(like='S2').values.flatten()
# LAI = field_data.filter(like='LAI').values.flatten()
CAB = field_data.filter(like='Cab').values.flatten() / 100.
# CAB = field_data.filter(like='VWC').values.flatten() /np.max(field_data.filter(like='VWC').values.flatten()) * (field_data.filter(like='Height').values.flatten() / 100.)**2

sm = field_data.filter(like='SM').values.flatten()
# sm = np.ones_like(LAI)*0.2

vv_field = field_data.filter(like='sigma_sentinel_vv').values.flatten()
vh_field = field_data.filter(like='sigma_sentinel_vh').values.flatten()
clay = 0.08
sand = 0.12
bulk = 1.5
f = 5.405

def run_SSRT(surface, vegetation, SAR):

    clay, sand, bulk, sm, s = surface
    d, LAI, coef, omega = vegetation
    f, theta = SAR

    ke = coef*np.sqrt(LAI)

    soil = Soil(mv=sm, s=s, clay=clay, sand=sand, f=f, bulk=bulk)

    can = OneLayer(canopy='turbid_isotropic', ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke)

    S = model.RTModel(surface=soil, canopy=can, models= {'surface': 'Oh92', 'canopy': 'turbid_isotropic'}, theta=theta, freq=f)
    S.sigma0()
    # pdb.set_trace()
    vv = S.__dict__['stot']['vv']
    vh = S.__dict__['stot']['hv']

    return vv, vh

def warper(s, m_vv, omega_vv, m_vh, omega_vh):

    coef_vv = m_vv * CAB
    coef_vh = m_vh * CAB
    SAR = [f, theta]

    surface_vv    = [clay, sand, bulk, sm, s]
    vegetation_vv = [1, LAI, coef_vv, omega_vv]


    surface_vh    = [clay, sand, bulk, sm, s]
    vegetation_vh = [1, LAI, coef_vh, omega_vh]

    vv = run_SSRT(surface_vv, vegetation_vv, SAR)[0]

    vh = run_SSRT(surface_vh, vegetation_vh, SAR)[1]

    return vv, vh

def cost(p):
    s , m_vv, omega_vv, m_vh, omega_vh = p
    #m_vv, m_vh = p[:len(LAI)], p[len(LAI):]
    #m_vv, m_vh = p
    vv, vh = warper(s, m_vv, omega_vv, m_vh, omega_vh)

    vv_diff = 10*np.log10(vv) - 10*np.log10(vv_field)
    vh_diff = 10*np.log10(vh) - 10*np.log10(vh_field)
    cost = vv_diff**2 + vh_diff**2

    if (~np.isfinite(cost)).all():
        cost = 999999
    else:
        mask = np.isfinite(cost)
        cost = np.nansum(cost[mask])
    print(cost)
    return cost

method = 'L-BFGS-B'

ps = 0.0105
pomega_vv = 0.027
pomega_vh = 0.0115
pm_vv = 1.4
pm_vh = 0.8
bounds = [(0.001, 10.03), (0.1,17), (0.001,10.07), (0.1,17), (0.001,10.04)]

guess = ps, pm_vv, pomega_vv, pm_vh, pomega_vh

res = minimize(cost, guess, method=method, bounds = bounds)

vv, vh = warper(res.x[0], res.x[1], res.x[2], res.x[3], res.x[4])
print(res.x)
plt.figure(figsize=(20,10))
plt.plot(10*np.log10(vv_field), label='vv dB in-situ')
plt.plot(10*np.log10(vh_field), label='vh dB in-situ')
plt.plot(10*np.log10(vv), label='vv dB calibrated model output')
plt.plot(10*np.log10(vh), label='vh dB calibrated model output')
plt.legend()
plt.ylabel('backscatter in dB')
plt.xlabel('time series')
plt.ylim((-27.5,-7.5))
plt.savefig('/media/tweiss/Daten/buffer/calibrationfieldbuffer50.png')

pdb.set_trace()

gp = gp_emulator.GaussianProcess(emulator_file='/media/tweiss/Daten/buffer/emulator_ssrt_s2_rfield.npz')

xxx = []
for i in range(len(vv)):
    xxx.append(random.uniform(0.5,1))




# prediction = np.column_stack((10*np.log10(vv_field), 10*np.log10(vh_field), LAI, CAB, np.rad2deg(theta)))
# prediction2 = np.column_stack((10*np.log10(vv), 10*np.log10(vh), LAI, CAB, np.rad2deg(theta)))
# prediction3 = np.column_stack((10*np.log10(vv)+xxx, 10*np.log10(vh)+xxx, LAI, CAB, np.rad2deg(theta)))
# prediction4 = np.column_stack((10*np.log10(vv)-xxx, 10*np.log10(vh)-xxx, LAI, CAB, np.rad2deg(theta)))

prediction = np.column_stack((10*np.log10(vv_field), 10*np.log10(vh_field),  np.rad2deg(theta)))
prediction2 = np.column_stack((10*np.log10(vv), 10*np.log10(vh),  np.rad2deg(theta)))
prediction3 = np.column_stack((10*np.log10(vv)+xxx, 10*np.log10(vh)+xxx,  np.rad2deg(theta)))
prediction4 = np.column_stack((10*np.log10(vv)-xxx, 10*np.log10(vh)-xxx,  np.rad2deg(theta)))



# prediction5 = np.column_stack((10*np.log10(vv)-1.5, 10*np.log10(vh)-1.5, LAI, CAB, np.rad2deg(theta)))

ypred, _, _ = gp.predict(prediction)
ypred2, _, _ = gp.predict(prediction2)
ypred3, _, _ = gp.predict(prediction3)
ypred4, _, _ = gp.predict(prediction4)
# ypred5, _, _ = gp.predict(prediction5)

plt.figure(figsize=(20,10))
plt.plot(sm, label='soil moisture (in-situ)')
plt.plot(ypred, label='soil moisture (retrieval with S1 backscatter)')
plt.plot(ypred2, label='soil moisture (retrieval with backscatter output of calibration)')
plt.plot(ypred3, '--', label='soil moisture (retrieval with backscatter output of calibration+random dB value between 0.5 and 1)')
plt.plot(ypred4, '--', label='soil moisture (retrieval with backscatter output of calibration+random dB value between -1 and -0.5)')
# plt.plot(ypred5, '--', label='soil moisture (retrieval with backscatter output of calibration+random dB value of 1)')
plt.ylim((0.1,0.4))
plt.title('Soil Moisture Retrieval Winter Wheat field 508')
plt.ylabel('soil moisture')
plt.xlabel('time series')
plt.legend()
plt.savefig('/media/tweiss/Daten/soil_moisture_retrieval2.png')
# plt.show()
plt.close()
plt.close()

sm = ypred
vv, vh = warper(res.x[0], res.x[1], res.x[2], res.x[3], res.x[4])
plt.plot(10*np.log10(vv))
plt.plot(10*np.log10(vh))
plt.plot(10*np.log10(vv_field))
plt.plot(10*np.log10(vh_field))

