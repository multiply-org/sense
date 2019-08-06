
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

def read_mni_data(path, file_name, extention, field, sep=';'):
    """ read MNI campaign data """
    df = pd.io.parsers.read_csv(os.path.join(path, file_name + extension), header=[0, 1], sep=sep)
    pd.to_datetime(df[field]['date']).dt.date
    df = df.set_index(pd.to_datetime(df[field]['date']).dt.date)
    # X=np.loadtxt("http://www2.geog.ucl.ac.uk/~ucfajlg/LMU_LAI_doy2017.txt")
    # xx = np.arange(1, 366)
    # lai_interp = np.interp(xx, X[:,0], X[:,1])
    # date_lai = datetime.date(2017, 1, 1) + datetime.timedelta(days=365)
    # xxx = datetime.date(2017,1,1) + xx * datetime.timedelta(days=1)
    # array = [[xxx],[lai_interp]]
    # tuples = list(zip(*arrays))
    # index = pd.MultiIndex.from_tuples(tuples, names=['first', 'second'])
    # dd = pd.DataFrame({'date':xxx, 'S2':lai_interp})
    # pdb.set_trace()
    dd = pd.io.parsers.read_csv(os.path.join('/media/tweiss/Daten/schlappi/params_S2.csv'), header=[0, 1], sep=',')
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

def smooth(x,window_len=11,window='hanning'):
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
    # pdb.set_trace()
    return df, df_agro, field_data, field_data_orbit, theta_field, sm_field, height_field, lai_field, vwc_field, pol_field, height_field2, lai_field2

### Optimization ###
#-----------------------------------------------------------------
def solve_fun(VALS):

    for i in range(len(var_opt)):
        dic[var_opt[i]] = VALS[i]

    ke = dic['coef'] * np.sqrt(dic['lai'])
    # ke = dic['coef'] * np.sqrt(dic['vwc'])
    # ke=1
    dic['ke'] = ke

    # surface
    soil = Soil(mv=dic['mv'], C_hh=dic['C_hh'], C_vv=dic['C_vv'], D_hh=dic['D_hh'], D_vv=dic['D_vv'], C_hv=dic['C_hv'], D_hv=dic['D_hv'], V2=dic['V2'], s=dic['s'], clay=dic['clay'], sand=dic['sand'], f=dic['f'], bulk=dic['bulk'], l=dic['l'])

    # canopy
    can = OneLayer(canopy=dic['canopy'], ke_h=dic['ke'], ke_v=dic['ke'], d=dic['d'], ks_h = dic['omega']*dic['ke'], ks_v = dic['omega']*dic['ke'], V1=dic['V1'], V2=dic['V2'], A_hh=dic['A_hh'], B_hh=dic['B_hh'], A_vv=dic['A_vv'], B_vv=dic['B_vv'], A_hv=dic['A_hv'], B_hv=dic['B_hv'])

    S = model.RTModel(surface=soil, canopy=can, models=models, theta=dic['theta'], freq=dic['f'])
    S.sigma0()

    return S.__dict__['stot'][pol[::-1]]

def fun_opt(VALS):
    # pdb.set_trace()

    # return(10.*np.log10(np.nansum(np.square(solve_fun(VALS)-dic['pol_value']))))
    return(np.nansum(np.square(solve_fun(VALS)-dic['pol_value'])))

def data_optimized_run(n, field_data, theta_field, sm_field, height_field, lai_field, vwc_field, pol):
    n = np.int(np.floor(n/2))

    if n > 0:
        field_data = field_data.drop(field_data.index[-n:])
        field_data = field_data.drop(field_data.index[0:n])
        theta_field = theta_field.drop(theta_field.index[-n:])
        theta_field = theta_field.drop(theta_field.index[0:n])

    sm_field = field_data.filter(like='SM')
    height_field = field_data.filter(like='Height')/100
    # height_field = field_data.filter(like='Cab')/100
    # lai_field = field_data.filter(like='S2')
    lai_field = field_data.filter(like='LAI')
    vwc_field = field_data.filter(like='VWC')

    vv_field = field_data.filter(like='sigma_sentinel_vv')
    vh_field = field_data.filter(like='sigma_sentinel_vh')

    pol_field = field_data.filter(like='sigma_sentinel_'+pol)
    return field_data, theta_field, sm_field, height_field, lai_field, vwc_field, vv_field, vh_field, pol_field
#-----------------------------------------------------------------

### Data preparation ###
#-----------------------------------------------------------------
# storage information
path = '/media/tweiss/Daten/schlappi'
file_name = 'multi10' # theta needs to be changed to for norm multi
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
bounds = [(0.001, 0.03), (0.1,7), (0.001,0.04), (0.1,7), (0.001,0.04)]

guess = ps, pm_vv, pomega_vv, pm_vh, pomega_vh
#guess = [1.4, 0.8]
res = minimize(cost, guess, method=method, bounds = bounds)
# solved = [0.02      , 0.00033854, 0.03      , 0.0105    , 0.01      ,
       # 0.0115    ]
# x = [0.38332326, 0.81561977, 0.45182801, 1.09407956, 1.23313943,
#        0.89714312, 1.17532486, 1.50653822, 2.74042195, 2.03820864,
#        1.2230573 , 0.60423899, 0.96138969, 0.91086824, 0.71054032,
#        0.77363155, 0.93829323, 0.73741354, 0.72120546, 0.76153766,
#        0.9678275 , 0.94895653, 0.9734852 , 0.99422879, 0.66862462,
#        0.77009955, 1.08004441, 0.66832872, 0.94710929, 1.19456028,
#        0.8068942 , 0.66108782, 1.16212088, 2.78515049, 0.91660068,
#        1.59904955, 0.80583734, 0.80024366, 2.57320055, 1.75382239,
#        0.9355827 , 0.96305127, 1.01815615, 0.72289155, 0.54808708,
#        1.0539223 , 0.817248  , 0.70664043, 1.10129883, 0.51504796,
#        0.61945984, 0.96577371, 0.62227635, 0.67872911, 1.77204458,
#        0.85231337, 0.66896275, 1.2463534 , 0.50000511, 0.06927394,
#        4.77100237, 3.91933208, 1.09205554, 1.00021858, 0.73680347,
#        1.02012947, 1.35574597, 0.66192406, 2.08511774, 1.64193998,
#        2.18341065, 1.98189336, 1.11971397, 0.55057736, 0.91173088,
#        0.7321862 , 0.54286355, 0.92850208, 0.81496359, 0.70728337,
#        0.79367924, 0.58072385, 0.47927211, 0.40612891, 0.61012018,
#        0.64781564, 0.6720015 , 0.6522005 , 0.52473976, 0.81640084,
#        0.45281958, 0.62125234, 0.67584209, 2.53014973, 1.14085623,
#        0.93231227, 2.48610563, 2.58956798, 0.54825612, 2.17920273,
#        0.86197714, 0.75560242, 0.58295085, 2.49279027, 0.82334317,
#        2.41240002, 1.0794502 , 1.32633502, 1.5532141 , 0.68328073,
#        0.70213215, 0.60101298, 0.36355202, 0.65086576, 0.41656313,
#        0.28684925, 1.58539464, 0.34905093, 0.1408389 , 2.05825474,
#        1.25328625, 1.11430827, 1.66874177, 8.68945532]
# xxx = cost([0.00548515, 0.73296915, 0.02264721, 0.00776594, 0.32751273,
       # 0.001     ])


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
plt.savefig('/media/tweiss/Daten/calibration.png')
# plt.show()
# pdb.set_trace()

gp = gp_emulator.GaussianProcess(emulator_file='/media/tweiss/Daten/emulator_ssrt_s2.npz')

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

