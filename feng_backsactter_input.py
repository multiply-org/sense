
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
    # field_data_orbit = filter_relativorbit(field_data, field, 95, 168)
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

theta_new = np.deg2rad(theta)

LAI_new = field_data.filter(like='S2').values.flatten()
LAI_new = field_data.filter(like='LAI').values.flatten()
CAB_new = field_data.filter(like='Cab').values.flatten() / 100.
# CAB = field_data.filter(like='VWC').values.flatten() /np.max(field_data.filter(like='VWC').values.flatten()) * (field_data.filter(like='Height').values.flatten() / 100.)**2

sm_field = field_data.filter(like='SM').values.flatten()
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

def wrapper(sm):

    s = 0.02706007
    m_vv = 0.68963347
    omega_vv = 0.02892507
    m_vh = 0.32105396
    omega_vh = 0.00370084

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
    sm = p
    #m_vv, m_vh = p[:len(LAI)], p[len(LAI):]
    #m_vv, m_vh = p
    vv, vh = wrapper(sm)

    vv_diff = 10*np.log10(vv) - 10*np.log10(vv_field_value)
    pdb.set_trace()
    # vh_diff = 10*np.log10(vh) - 10*np.log10(vh_field_value)
    cost = vv_diff**2

    if (~np.isfinite(cost)).all():
        cost = 999999
    else:
        mask = np.isfinite(cost)
        cost = np.nansum(cost[mask])
    print(cost)
    return cost

method = 'L-BFGS-B'

psm = 0.2
bounds = [(0.001, 0.5)]

soilmoisture = []

# for i, vv_field_value in enumerate(vv_field):
#     LAI = LAI_new[i]
#     CAB = CAB_new[i]
#     theta = theta_new[i]
#     guess = psm
#     res = minimize(cost, guess, method=method, bounds = bounds)
#     soilmoisture.append(res.x[0])

lookup = {}
lookup_lai = {}
lookup_cab = {}
lookup_theta = {}
lookup_vv = {}
vv_dreck = []


# index =

columns = ['LAI', 'CAB', 'theta', 'vv', 'sm']

dx = pd.DataFrame(columns=columns)
j=0
for i, LAI in enumerate(np.arange(0.1,7,0.01)):
    for ii, CAB in enumerate(np.arange(0.1,1,0.01)):
        for iii, theta in enumerate(np.arange(0.1, 1, 0.01)):
            for iiii, sm in enumerate(np.arange(0.01, 0.5, 0.01)):
                vv, vh = wrapper(sm)

                vv_dreck.append(vv)
                dx2 = pd.DataFrame([[LAI,CAB,theta,float(vv),sm]], columns=columns,index=[j])
                dx =dx.append(dx2)
                j =j+1

                # lookup[str(vv),str(LAI),str(CAB),str(theta)] = sm
                # lookup_lai[str(LAI)] = LAI
                # lookup_cab[str(CAB)] = CAB
                # lookup_theta[str(theta)] = theta
                # lookup_vv[str(vv)] = vv

dx.to_csv('/media/tweiss/Daten/schlappi/lookup.csv')
pdb.set_trace()

index_LAI = abs(dx['LAI'] - 0.2).idxmin()
dx_new = dx[dx['LAI']==dx['LAI'][index_LAI]]

index_CAB = abs(dx_new['CAB'] - 0.2).idxmin()
dx_new = dx_new[dx_new['CAB']==dx_new['CAB'][index_CAB]]

index_theta = abs(dx_new['theta'] - 0.2).idxmin()
dx_new = dx_new[dx_new['theta']==dx_new['theta'][index_theta]]

index_vv = abs(dx_new['vv'] - 0.2).idxmin()
dx_new = dx_new[dx_new['vv']==dx_new['vv'][index_vv]]

pdb.set_trace()

# xv, yv, zv, vx, se, sv = np.meshgrid(np.arange(0.1,7,1), np.arange(0.1,1,0.1), np.arange(0.1, 1, 0.1), vv_dreck, np.arange(0.01, 0.5, 0.1), 1)
# np.min((abs(abs(dx['LAI'])-0.2)+0.2

for i, LAI in enumerate(np.arange(0.1,7,1)):
    for ii, CAB in enumerate(np.arange(0.1,1,0.1)):
        for iii, theta in enumerate(np.arange(0.1, 1, 0.1)):
            for iiii, sm in enumerate(np.arange(0.01, 0.5, 0.1)):
                vv, vh = wrapper(sm)
                coord=[LAI,CAB,theta,float(vv),sm]
                index = np.argwhere((xv==coord[0]) & (yv==coord[1]) & (zv==coord[2]) & (vx==coord[3]) & (se==coord[4]))
                print(index)
                sv[index] = sm




coord=[0.1,0.1,0.1,0.10732377322886114]
np.argwhere((xv==coord[0]) & (yv==coord[1]) & (zv==coord[2]) & (vx==coord[3]))
