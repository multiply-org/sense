
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
    dd = pd.io.parsers.read_csv(os.path.join('/media/wodan/Volume/Arbeit/params_S2.csv'), header=[0, 1], sep=',')
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
    height_field2 = field_data.filter(like='Height')/100
    height_field = field_data.filter(like='Cab')/100
    lai_field2 = field_data.filter(like='LAI')
    lai_field = field_data.filter(like='S2')
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
    # height_field = field_data.filter(like='Height')/100
    height_field = field_data.filter(like='Cab')/100
    lai_field = field_data.filter(like='S2')
    vwc_field = field_data.filter(like='VWC')

    vv_field = field_data.filter(like='sigma_sentinel_vv')
    vh_field = field_data.filter(like='sigma_sentinel_vh')

    pol_field = field_data.filter(like='sigma_sentinel_'+pol)
    return field_data, theta_field, sm_field, height_field, lai_field, vwc_field, vv_field, vh_field, pol_field
#-----------------------------------------------------------------

### Data preparation ###
#-----------------------------------------------------------------
# storage information
path = '/media/wodan/Volume/Arbeit'
file_name = 'multi_10_neu' # theta needs to be changed to for norm multi
extension = '.csv'

path_agro = '/media/nas_data/2017_MNI_campaign/field_data/meteodata/agrarmeteorological_station'
file_name_agro = 'Eichenried_01012017_31122017_hourly'
extension_agro = '.csv'

field = '508_high'
pol = 'vv'

df, df_agro, field_data, field_data_orbit, theta_field, sm_field, height_field, lai_field, vwc_field, pol_field, height_field2, lai_field2 = read_data(path, file_name, extension, field, path_agro, file_name_agro, extension_agro)


# X=np.loadtxt("http://www2.geog.ucl.ac.uk/~ucfajlg/LMU_Kaska_params.txt", delimiter=',')
# xx = np.arange(1, 366)
# lai_interp = np.interp(xx, X[:,0], X[:,1])
# Cab_interp = np.interp(xx, X[:,0], X[:,2])
# Cbrown_interp = np.interp(xx, X[:,0], X[:,3])


# date_lai = datetime.date(2017, 1, 1) + datetime.timedelta(days=365)
# xxx = datetime.date(2017,1,1) + xx * datetime.timedelta(days=1)

# dd = pd.DataFrame({'date':xxx, 'LAI':lai_interp, 'Cab':Cab_interp, 'Cbrown':Cbrown_interp})
# pdb.set_trace()




#-----------------------------------------------------------------

### Settings SenSe module ###
#-----------------------------------------------------------------
## Choose models
#---------------
surface = 'Oh92'
# surface = 'Oh04'
# surface = 'Dubois95'
# surface = 'WaterCloud'
# surface = 'I2EM'
# canopy = 'turbid_isotropic'
# canopy = 'turbid_rayleigh'
canopy = 'water_cloud'

models = {'surface': surface, 'canopy': canopy}

## Parameters
#----------------
# anyway
freq = 5.405

### Surface
# Oh92 (freq, clay, sand, bulk, s, theta, mv)
# Oh04 (freq, s, theta, mv)
# Dubois95 (freq, clay, sand, bulk, s, theta, mv)
# Water Cloud (C, D)
#-------------------------------------------------
clay = 0.08
sand = 0.12
bulk = 1.5
# s = 0.0105
s = 0.007

C_hh = 0
D_hh = 0
C_hv = -22.5
D_hv = 4
C_vv = -14.609339
D_vv = 12.884086

### Canopy
# Water Cloud (A, B, V1, V2, theta)
# SSRT (coef, omega, theta)
#-----------------------------------
A_hh = 0
B_hh = 0
A_hv = 0.029
B_hv = 0.0013
A_vv = 0.0029
B_vv = 0.13
V1 = lai_field.values.flatten()
V2 = V1 # initialize in surface model
coef = 1.
omega = 0.027
omega = 0.115
# IEM
l = 0.01
#-----------------------------------------------------------------

### Run SenSe module
#-----------------------------------------------------------------
#### Choose models
#-----------------
# surface = 'Oh92'
# surface = 'Oh04'
surface = 'Dubois95'
# surface = 'WaterCloud'
surface = 'I2EM'
canopy = 'turbid_isotropic'
# canopy = 'turbid_rayleigh'
# canopy = 'water_cloud'

models = {'surface': surface, 'canopy': canopy}

opt_mod = 'time invariant'
opt_mod = 'time variant'

surface_list = ['Oh92']
canopy_list = ['turbid_isotropic']

fig, ax = plt.subplots(figsize=(20, 10))
# plt.title('Winter Wheat')
plt.ylabel('Backscatter [dB]', fontsize=15)
plt.xlabel('Date', fontsize=15)
plt.tick_params(labelsize=12)
ax.set_ylim([-22.5,-7.5])

colormaps = ['Greens', 'Purples', 'Blues', 'Oranges', 'Reds', 'Greys', 'pink', 'bone']
j = 0

colormap = plt.get_cmap(colormaps[j])
colors = [colormap(jj) for jj in np.linspace(0.35, 1., 3)]

for k in surface_list:

    for kk in canopy_list:
        df, df_agro, field_data, field_data_orbit, theta_field, sm_field, height_field, lai_field, vwc_field, pol_field, height_field2, lai_field2 = read_data(path, file_name, extension, field, path_agro, file_name_agro, extension_agro)
        clay = 0.08
        sand = 0.12
        bulk = 1.5
        # s = 0.0105
        s = 0.007

        C_hh = 0
        D_hh = 0
        C_hv = -22.5
        D_hv = 4
        C_vv = -14.609339
        D_vv = 12.884086

        ### Canopy
        # Water Cloud (A, B, V1, V2, theta)
        # SSRT (coef, omega, theta)
        #-----------------------------------
        A_hh = 0
        B_hh = 0
        A_hv = 0.029
        B_hv = 0.0013
        A_vv = 0.0029
        B_vv = 0.13
        V1 = lai_field.values.flatten()
        V2 = V1 # initialize in surface model
        coef = 1.
        omega = 0.027
        # omega = 0.0115
        # IEM
        l = 0.01


        surface = k
        canopy = kk
        models = {'surface': surface, 'canopy': canopy}

        #### Optimization
        #-----------------

        if opt_mod == 'time invariant':

            dic = {"mv":sm_field.values.flatten(), "C_hh":C_hh, "C_vv":C_vv, "D_hh":D_hh, "D_vv":D_vv, "C_hv":C_hv, "D_hv":D_hv, "s":s, "clay":clay, "sand":sand, "f":freq, "bulk":bulk, "l":l, "canopy":canopy, "d":height_field.values.flatten(), "V1":V1, "V2":V2, "A_hh":A_hh, "B_hh":B_hh, "A_vv":A_vv, "B_vv":B_vv, "A_hv":A_hv, "B_hv":B_hv, "lai":lai_field.values.flatten(), "vwc":vwc_field.values.flatten(), "pol_value":pol_field.values.flatten(), "theta":theta_field.values.flatten(), "omega": omega, "coef": coef}

            if canopy == 'turbid_isotropic':
                var_opt = ['coef']
                guess = [2.]
                bounds = [(0.,5.5)]
            elif surface == 'WaterCloud' and canopy == 'water_cloud':
                var_opt = ['A_vv', 'B_vv', 'A_hv', 'B_hv', 'C_vv', 'D_vv', 'C_hv', 'D_hv']
                guess = [A_vv, B_vv, A_hv, B_hv, C_vv, D_vv, C_hv, D_hv]
                bounds = [(0.,1), (0.,1), (0.,1), (0.,1), (-20.,-1.), (1.,20.), (-20.,-1.), (1.,20.)]
            elif canopy == 'water_cloud':
                var_opt = ['A_vv', 'B_vv', 'A_hv', 'B_hv']
                guess = [A_vv, B_vv, A_hv, B_hv]
                bounds = [(0.,1), (0.,1), (0.,1), (0.,1)]

            method = 'L-BFGS-B'

            res = minimize(fun_opt,guess,bounds=bounds, method=method)

            fun_opt(res.x)
            aaa = res.x

        if opt_mod == 'time variant':
            aaa = [[],[],[],[],[],[],[],[],[],[],[],[]]
            n=7

            for i in range(len(pol_field.values.flatten())-n+1):

                if type(coef) == float:
                    dic = {"mv":sm_field.values.flatten()[i:i+n], "C_hh":C_hh, "C_vv":C_vv, "D_hh":D_hh, "D_vv":D_vv, "C_hv":C_hv, "D_hv":D_hv, "V2":V2[i:i+n], "s":s, "clay":clay, "sand":sand, "f":freq, "bulk":bulk, "l":l, "canopy":canopy, "d":height_field.values.flatten()[i:i+n], "V1":V1[i:i+n], "A_hh":A_hh, "B_hh":B_hh, "A_vv":A_vv, "B_vv":B_vv, "A_hv":A_hv, "B_hv":B_hv, "lai":lai_field.values.flatten()[i:i+n], "vwc":vwc_field.values.flatten()[i:i+n], "pol_value":pol_field.values.flatten()[i:i+n], "theta":theta_field.values.flatten()[i:i+n], "omega": omega, "coef": coef}
                else:
                    dic = {"mv":sm_field.values.flatten()[i:i+n], "C_hh":C_hh, "C_vv":C_vv, "D_hh":D_hh, "D_vv":D_vv, "C_hv":C_hv, "D_hv":D_hv, "V2":V2[i:i+n], "s":s, "clay":clay, "sand":sand, "f":freq, "bulk":bulk, "l":l, "canopy":canopy, "d":height_field.values.flatten()[i:i+n], "V1":V1[i:i+n], "A_hh":A_hh, "B_hh":B_hh, "A_vv":A_vv, "B_vv":B_vv, "A_hv":A_hv, "B_hv":B_hv, "lai":lai_field.values.flatten()[i:i+n], "vwc":vwc_field.values.flatten()[i:i+n], "pol_value":pol_field.values.flatten()[i:i+n], "theta":theta_field.values.flatten()[i:i+n], "omega": omega, "coef": coef[i:i+n]}

                if canopy == 'turbid_isotropic' and surface == 'WaterCloud':
                    var_opt = ['coef', 'C_vv', 'D_vv', 'C_hv', 'D_hv']
                    guess = [1., C_vv, D_vv, C_hv, D_hv]
                    bounds = [(0.1,5.5), (-20.,-1.), (1.,20.), (-20.,-1.), (1.,20.)]
                elif canopy == 'turbid_isotropic':
                    var_opt = ['coef']
                    guess = [0.1]
                    bounds = [(0.,5)]
                elif surface == 'WaterCloud' and canopy == 'water_cloud':
                    var_opt = ['A_vv', 'B_vv', 'A_hv', 'B_hv', 'C_vv', 'D_vv', 'C_hv', 'D_hv']
                    guess = [A_vv, B_vv, A_hv, B_hv, C_vv, D_vv, C_hv, D_hv]
                    bounds = [(0.,5), (guess[1]*0.55, guess[1]*1.55), (0.,5), (guess[3]*0.75, guess[3]*1.25), (-20.,-1.), (1.,20.), (-20.,-1.), (1.,20.)]
                    # var_opt = ['C_vv', 'D_vv', 'C_hv', 'D_hv']
                    # guess = [C_vv, D_vv, C_hv, D_hv]
                    # bounds = [(-20.,-1.), (1.,20.), (-20.,-1.), (1.,20.)]
                elif canopy == 'water_cloud':
                    var_opt = ['A_vv', 'B_vv', 'A_hv', 'B_hv']
                    guess = [A_vv, B_vv, A_hv, B_hv]
                    bounds = [(0.,1), (0.,1), (0.,1), (0.,1)]

                # var_opt = ['omega']
                # guess = [0.1]
                # bounds = [(0.,5.5)]

                # var_opt = ['s', 'coef', 'omega']
                # guess = [0.01, 0.1, 0.01]
                # bounds = [(0.001,0.03),(0.,2.5),(0.001,0.1)]

                # var_opt = ['C_hv', 'D_hv']
                # guess = [-13, 14]
                # bounds = [(-200.,100.),(-200.,400.)]

                # var_opt = ['A_vv', 'B_vv']

                # try:
                #     guess = [res.x[0], res.x[1]]
                # except:
                #     guess = [0.005, 0.09]
                # # bounds = [(0.000,5.),(0.001,5.)]
                # bounds = [(guess[0]*0.75, guess[0]*1.25), (guess[1]*0.75, guess[1]*1.25)]
                # bounds = [(guess[0]*0.9, guess[0]*1.1), (guess[1]*0.75, guess[1]*1.25)]
                # var_opt = ['coef', 'omega']
                # guess = [0.1, 0.22]
                # bounds = [(0.,5.5),(0.00001,0.2)]
                method = 'L-BFGS-B'
                # method = 'trust-exact'

                res = minimize(fun_opt,guess,bounds=bounds, method=method)

                fun_opt(res.x)

                for j in range(len(res.x)):
                    aaa[j].append(res.x[j])

            field_data, theta_field, sm_field, height_field, lai_field, vwc_field, vv_field, vh_field, pol_field = data_optimized_run(n, field_data, theta_field, sm_field, height_field, lai_field, vwc_field, pol)
            V1 = lai_field.values.flatten()
            V2 = V1 # initialize in surface model
            # pdb.set_trace()

        #-----------------------------------------------------------------

        for i in range(len(res.x)):
            exec('%s = %s' % (var_opt[i],aaa[i]))

        ke = coef * np.sqrt(lai_field.values.flatten())
        # ke = smooth(ke, 11)

        soil = Soil(mv=sm_field.values.flatten(), C_hh=np.array(C_hh), C_vv=np.array(C_vv), D_hh=np.array(D_hh), D_vv=np.array(D_vv), C_hv=np.array(C_hv), D_hv=np.array(D_hv), s=s, clay=clay, sand=sand, f=freq, bulk=bulk, l=l)

        can = OneLayer(canopy=canopy, ke_h=ke, ke_v=ke, d=height_field.values.flatten(), ks_h = omega*ke, ks_v = omega*ke, V1=np.array(V1), V2=np.array(V2), A_hh=np.array(A_hh), B_hh=np.array(B_hh), A_vv=np.array(A_vv), B_vv=np.array(B_vv), A_hv=np.array(A_hv), B_hv=np.array(B_hv))

        S = model.RTModel(surface=soil, canopy=can, models=models, theta=theta_field.values.flatten(), freq=freq)
        S.sigma0()
#-----------------------------------------------------------------
        date = field_data.index

        colormap = plt.get_cmap(colormaps[j])
        colors = [colormap(jj) for jj in np.linspace(0.35, 1., 4)]

        # ax.plot(10*np.log10(pol_field), 'ks-', label='Sentinel-1 Pol: ' + pol, linewidth=3)
        # ax.plot(date, 10*np.log10(S.__dict__['s0g'][pol[::-1]]), color=colors[0], marker='s', linestyle='--', label=pol+' s0g')
        # ax.plot(date, 10*np.log10(S.__dict__['s0c'][pol[::-1]]), color=colors[1], marker='s', linestyle='--', label=pol+' s0c')
        # ax.plot(date, 10*np.log10(S.__dict__['s0cgt'][pol[::-1]]), 'ms-', label=pol+' s0cgt')
        # ax.plot(date, 10*np.log10(S.__dict__['s0gcg'][pol[::-1]]), 'ys-', label=pol+' s0gcg')

        mask = ~np.isnan(pol_field.values.flatten()) & ~np.isnan(S.__dict__['stot'][pol[::-1]])
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress((pol_field.values.flatten()[mask]), (S.__dict__['stot'][pol[::-1]][mask]))
        slope1, intercept1, r_value1, p_value1, std_err1 = scipy.stats.linregress(10*np.log10(pol_field.values.flatten())[mask], 10*np.log10(S.__dict__['stot'][pol[::-1]])[mask])
        rmse = rmse_prediction(10*np.log10(pol_field.values.flatten()), 10*np.log10(S.__dict__['stot'][pol[::-1]]))

        if k == 'Oh92':
            colors = 'blue'
        elif k == 'Oh04':
            colors = 'red'
        elif k == 'Dubois95':
            colors = 'orange'
        elif k == 'WaterCloud':
            colors = 'purple'


        if kk == 'turbid_isotropic':

            ax.plot(date, 10*np.log10(S.__dict__['stot'][pol[::-1]]), color=colors, marker='s', linestyle='dashed', label=S.models['surface']+ ' + ' +  S.models['canopy'] + ' Pol: ' + pol + '; RMSE: ' + str(rmse)[0:4] + '; R2: ' + str(r_value)[0:4] + ' ' + str(r_value1)[0:4])
        else:
            ax.plot(date, 10*np.log10(S.__dict__['stot'][pol[::-1]]), color=colors, marker='s', label=S.models['surface']+ ' + ' +  S.models['canopy'] + ' Pol: ' + pol + '; RMSE: ' + str(rmse)[0:4] + '; R2: ' + str(r_value)[0:4] + ' ' + str(r_value1)[0:4])

        j = j+1
ax.plot(10*np.log10(pol_field), 'ks-', label='Sentinel-1 Pol: ' + pol, linewidth=3)
plt.legend()
plt.title('Winter Wheat 508')

ax1 = ax.twinx()

B = height_field.values.flatten()*ke

ax1.plot(pd.to_datetime(date), B*2, label='B')
ax1.plot(pd.to_datetime(date), np.sqrt(lai_field.values.flatten()), label='LAI')
ax1.plot(pd.to_datetime(date), coef, label='coef')
ax1.plot(pd.to_datetime(date), ke, label='ke')
ax1.plot(pd.to_datetime(date), height_field.values.flatten(), label='height')
# ax1.plot(pd.to_datetime(date), A_vv, label='A')
# ax1.plot(pd.to_datetime(date), B_vv, label='B')
# ax1.plot(pd.to_datetime(date), C_vv, label='C')
# ax1.plot(pd.to_datetime(date), D_vv, label='D')
# ax1.plot(pd.to_datetime(dd.date),lai_interp*Cab_interp/20, label='lai*cab')

# ax1.plot(pd.to_datetime(dd.date),lai_interp, label='lai')
# ax1.plot(pd.to_datetime(dd.date),Cab_interp/10., label='cab')
# ax1.plot(pd.to_datetime(dd.date),Cbrown_interp*10, label='cbrown')
ax1.set_xlim(['2017-03-20', '2017-07-20'])
ax.set_ylim([-25,-10])
plt.legend()
plt.savefig('/media/wodan/Volume/Arbeit/vv_'+opt_mod)
pdb.set_trace()

































































colors = ['gs-', 'rs-', ]

##### Plot

date = field_data.index

# fig, ax = plt.subplots(figsize=(20, 10))
# # plt.title('Winter Wheat')
# plt.ylabel('Backscatter [dB]', fontsize=15)
# plt.tick_params(labelsize=12)

ax.plot(10*np.log10(pol_field), 'ks-', label='Sentinel-1 Pol: ' + pol, linewidth=3)
ax.plot(date, 10*np.log10(S.__dict__['s0g'][pol[::-1]]), 'rs-', label=pol+' s0g')
ax.plot(date, 10*np.log10(S.__dict__['s0c'][pol[::-1]]), 'cs-', label=pol+' s0c')
# ax.plot(date, 10*np.log10(S.__dict__['s0cgt'][pol[::-1]]), 'ms-', label=pol+' s0cgt')
# ax.plot(date, 10*np.log10(S.__dict__['s0gcg'][pol[::-1]]), 'ys-', label=pol+' s0gcg')
ax.plot(date, 10*np.log10(S.__dict__['stot'][pol[::-1]]), 'C1s-', label=S.models['surface']+ ' + ' +  S.models['canopy'] + ' Pol: ' + pol)
ax.legend()
ax.legend(loc=2, fontsize=12)
ax.set_ylim([-35,-5])
ax.set_xlim(['2017-03-20', '2017-07-20'])

mask = ~np.isnan(pol_field.values.flatten()) & ~np.isnan(S.__dict__['stot'][pol[::-1]])
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress((pol_field.values.flatten()[mask]), (S.__dict__['stot'][pol[::-1]][mask]))
slope1, intercept1, r_value1, p_value1, std_err1 = scipy.stats.linregress(10*np.log10(pol_field.values.flatten())[mask], 10*np.log10(S.__dict__['stot'][pol[::-1]])[mask])
rmse = rmse_prediction(10*np.log10(pol_field.values.flatten()), 10*np.log10(S.__dict__['stot'][pol[::-1]]))

plt.title('Winter Wheat, R2 = ' + str(r_value) + ' RMSE = ' + str(rmse))

pdb.set_trace()


































#-----------------------------------------------------------------


### Run SenSe module
#-----------------------------------------------------------------

V1 = lai_field.values.flatten()
V2 = lai_field.values.flatten()

soil = Soil(mv=sm_field.values.flatten(), C_hh=C_hh, C_vv=C_vv, D_hh=D_hh, D_vv=D_vv, C_hv=C_hv, D_hv=D_hv, V2=V2, s=s, clay=clay, sand=sand, f=freq, bulk=bulk, l=l)


can = OneLayer(canopy=canopy, ke_h=ke, ke_v=ke, d=height_field.values.flatten(), ks_h = omega*ke, ks_v = omega*ke, V1=V1, V2=V2, A_hh=A_hh, B_hh=B_hh, A_vv=A_vv, B_vv=B_vv, A_hv=A_hv, B_hv=B_hv)

S = model.RTModel(surface=soil, canopy=can, models=models, theta=theta_field.values.flatten(), freq=freq)
S.sigma0()

if field == '508_med':
    S.__dict__['stot'][pol[::-1]][13:16] = np.nan


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
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==44]), 'ys-', label=44)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==117]), 'ms-', label=117)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==95]), 'rs-', label=95)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==168]), 'gs-', label=168)

ax.plot(date, 10*np.log10(S.__dict__['s0g'][pol[::-1]]), 'rs-', label=pol+' s0g')
ax.plot(date, 10*np.log10(S.__dict__['s0c'][pol[::-1]]), 'cs-', label=pol+' s0c')
# ax.plot(date, 10*np.log10(S.__dict__['s0cgt'][pol[::-1]]), 'ms-', label=pol+' s0cgt')
# ax.plot(date, 10*np.log10(S.__dict__['s0gcg'][pol[::-1]]), 'ys-', label=pol+' s0gcg')
ax.plot(date, 10*np.log10(S.__dict__['stot'][pol[::-1]]), 'C1s-', label=S.models['surface']+ ' + ' +  S.models['canopy'] + ' Pol: ' + pol)
ax.legend()
ax.legend(loc=2, fontsize=12)

# ax2 = ax.twinx()
# ax2.plot(date, aaa, label='fitted parameter')
# ax2.plot(lai_field, color='green', label='LAI')
# ax2.plot(vwc_field, color='blue', label='VWC')
# ax2.legend(loc=1, fontsize=12)
# ax2.tick_params(labelsize=12)
# ax2.set_ylabel('LAI [m$^3$/m$^3$]', fontsize=15, color='green')
ax3 = ax.twinx()
# ax3.set_ylabel('VWC [kg/m$^2$]', fontsize=15, color='blue')
ax3.tick_params(labelsize=12)
# lns2 = ax3.plot(date,ke, color='blue', label='Volume extinction coefficient')
# ax3.set_ylabel('Volume extinction coefficient [Np/m]', fontsize=15)
# ax3.yaxis.label.set_color('blue')
# ax3.plot(date,ke, 'bs-')
# ax3.plot(date,coef, 'ys-')

# ax3.plot(date[field_data[(field,'relativeorbit')]==44], ke[field_data[(field,'relativeorbit')]==44], 'ys', label=44)
# ax3.plot(date[field_data[(field,'relativeorbit')]==117], ke[field_data[(field,'relativeorbit')]==117], 'ms', label=117)
# ax3.plot(date[field_data[(field,'relativeorbit')]==95], ke[field_data[(field,'relativeorbit')]==95], 'rs', label=95)
# ax3.plot(date[field_data[(field,'relativeorbit')]==168], ke[field_data[(field,'relativeorbit')]==168], 'gs', label=168)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==117]), 'ms-', label=117)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==95]), 'rs-', label=95)
# ax.plot(10*np.log10(pol_field[field_data[(field,'relativeorbit')]==168]), 'gs-', label=168)



# ax3.plot(date,omega, 'ms-')
# ax4 = ax2.twinx()
# ax4.plot(sm_field)
# ax4.set_ylim([-0.8,0.4])

# ax5 = ax2.twinx()
# ax5.plot(height_field, 'g-')

ax6 = ax.twinx()
ax6.tick_params(labelsize=12)
lns1 = ax6.plot(lai_field, 'g', label='VWC')
ax6.set_ylabel('LAI [/]', fontsize=15, color='green')
ax6.yaxis.label.set_color('green')
ax6.spines['right'].set_position(('outward', 60))

# lns = lns1+lns2
# labs = [l.get_label() for l in lns]
# ax3.legend(lns, labs, loc=0, fontsize=12)



# ax7 = ax2.twinx()
# ax7.plot(lai_field, 'c')

# ax8 = ax.twinx()
# ax8.plot(sm_field, '+')
# ax8.set_ylabel('SM')
# ax8.spines['right'].set_position(('outward', 30))


# ax9 = ax.twinx()
# ax9.bar(df_agrar.index, df_agrar['SUM_NN050'], 0.1)
# ax9.set_ylabel('precip', ha='left')
# ax9.spines['right'].set_position(('outward', 60))
# ax9.set_ylim([0,15])

ax.grid(linestyle='-', linewidth=1)
ax.grid(b=True, which='minor', linestyle='--', linewidth=0.5)

days = mdates.DayLocator()
ax.xaxis.set_minor_locator(days)

months = MonthLocator()
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b %Y"))

# ax.set_ylim([np.min(10*np.log10(pol_field.values.flatten()))-2, np.max(10*np.log10(pol_field.values.flatten()))+2])
ax.set_ylim([-30,-5])
# ax2.get_yaxis().set_ticks([])
# ax2.set_ylim([0,8])
# ax3.set_ylim([0,8])
ax.set_xlim(['2017-03-20', '2017-07-15'])
# ax.set_xlim(['2017-06-01', '2017-09-30'])

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress((pol_field.values.flatten()), (S.__dict__['stot'][pol[::-1]]))
slope1, intercept1, r_value1, p_value1, std_err1 = scipy.stats.linregress(10*np.log10(pol_field.values.flatten()), 10*np.log10(S.__dict__['stot'][pol[::-1]]))
rmse = rmse_prediction(10*np.log10(pol_field.values.flatten()), 10*np.log10(S.__dict__['stot'][pol[::-1]]))

plt.title('Winter Wheat, R2 = ' + str(r_value) + ' RMSE = ' + str(rmse))
plt.savefig('/media/tweiss/Daten/plots/plot_'+field+'_'+pol+'_'+file_name+'_'+S.models['surface']+'_'+S.models['canopy'])
plt.close()


pdb.set_trace()

# ke = np.array([ 2.06326886,  1.91916714,  1.77367103,  1.63778107,  1.528702  ,
#         1.46604885,  1.46253485,  1.51355678,  1.60051594,  1.69869197,
#         1.78462896,  1.84024706,  1.85075172,  1.81189095,  1.73825148,
#         1.65175839,  1.57682432,  1.53837552,  1.55148688,  1.62179581,
#         1.74169402,  1.88590305,  2.02932635,  2.14663887,  2.21078958,
#         2.21276201,  2.15754037,  2.06089996,  1.94999007,  1.84663998,
#         1.76203146,  1.69526738,  1.6298006 ,  1.55229176,  1.46191104,
#         1.36124346,  1.25958148,  1.16926533,  1.09563444,  1.03437644,
#         0.97349969,  0.90241684,  0.81561377,  0.71172699,  0.59756655,
#         0.48556321,  0.3871495 ,  0.30900663,  0.25081152,  0.20861753,
#         0.17712126,  0.14943961,  0.12191173,  0.09462769,  0.06671757,
#         0.03757646])
# C_vv = np.array([-14.62021844, -14.52733754,  -8.07220241, -14.74143993,
#        -18.42041188, -14.36304957, -14.28024849, -14.3342842 ,
#        -14.36442823, -14.23628017,  -9.25625934,  -8.19237892,
#         -5.28158493,  -3.87026981,  -3.45584735,  -3.98370261,
#        -13.18425257, -12.98646127, -12.92703971, -13.04234691,
#        -11.15906684, -11.15784566, -11.21186719, -11.5211018 ,
#        -11.308834  , -11.19225841, -11.31916709, -10.84673287,
#        -10.91784905, -12.99996228, -12.99996436, -12.99997557,
#        -12.99998938, -13.        , -12.99998988, -12.99998993,
#        -13.        , -13.        , -12.99998724, -12.99998405,
#        -12.99997974, -12.99997786, -12.99997295, -12.99996202,
#         -3.49934522,  -3.92486052,  -3.58896055,  -3.51956506,
#         -4.75782816,  -3.60160112,  -2.99269511,   0.        ,
#          0.        ,  -0.66617689,   0.        ,   0.        ])
# D_vv = np.array([ 13.13525467,  13.19322237,   1.37567403,  13.0790507 ,
#         20.        ,  13.26559902,  13.2763212 ,  13.22528786,
#         13.18792346,  13.20142913,   6.05404625,   4.56398804,
#          0.15512926,   0.        ,   0.        ,   1.21771858,
#         20.        ,  20.        ,  20.        ,  20.        ,
#         15.00725676,  14.96739904,  14.93159261,  14.77349734,
#         14.88187912,  14.94077546,  14.86500294,  15.02409645,
#         14.98689808,  14.00002015,  14.00001897,  14.00001251,
#         14.00000518,  14.        ,  14.00000513,  14.00000512,
#         14.        ,  14.        ,  14.00000639,  14.0000078 ,
#         14.0000097 ,  14.00000997,  14.0000117 ,  14.00001584,
#         17.71884858,  17.59050298,  20.        ,  20.        ,
#         20.        ,  18.39475902,  18.06682578,  12.51381333,
#         13.69266562,  15.42614678,  11.33772659,  11.57183062])
aaa = []
bbb = []
ccc = []

field_plot = ['508_high', '508_low', '508_med']
i = 0
colors = ['ks', 'ys', 'ms', 'rs']

for field in field_plot:
    df, df_agro, field_data, field_data_orbit, theta_field, sm_field, height_field, lai_field, vwc_field, pol_field = read_data(path, file_name, extension, field, path_agro, file_name_agro, extension_agro)
    field_data, theta_field, sm_field, height_field, lai_field, vwc_field, vv_field, vh_field, pol_field = data_optimized_run(n, field_data, theta_field, sm_field, height_field, lai_field, vwc_field, pol)

    soil = Soil(mv=sm_field.values.flatten(), C_hh=C_hh, C_vv=C_vv, D_hh=D_hh, D_vv=D_vv, C_hv=C_hv, D_hv=D_hv, V2=lai_field.values.flatten(), s=s, clay=clay, sand=sand, f=freq, bulk=bulk, l=l)


    can = OneLayer(canopy=canopy, ke_h=ke, ke_v=ke, d=height_field.values.flatten(), ks_h = omega*ke, ks_v = omega*ke, V1=lai_field.values.flatten(), V2=lai_field.values.flatten(), A_hh=A_hh, B_hh=B_hh, A_vv=A_vv, B_vv=B_vv, A_hv=A_hv, B_hv=B_hv)

    S = model.RTModel(surface=soil, canopy=can, models=models, theta=theta_field.values.flatten(), freq=freq)
    S.sigma0()

    plt.plot(10*np.log10(pol_field.values.flatten()),10*np.log10(S.__dict__['stot'][pol[::-1]]), colors[i], label=field)

    aaa = np.append(aaa, 10*np.log10(pol_field.values.flatten()))
    bbb = np.append(bbb, 10*np.log10(S.__dict__['stot'][pol[::-1]]))
    i = i+1


plt.ylabel(surface + ' ' + canopy + ' [dB]')
plt.xlabel('Sentinel-1 [dB]')
plt.legend()
x = np.linspace(np.min(10*np.log10(pol_field.values.flatten()))-2, np.max(10*np.log10(pol_field.values.flatten()))+2, 16)
plt.plot(x,x)
www = rmse_prediction(aaa, bbb)
# slope, intercept, r_value, p_value, std_err = linregress(aaa[~np.isnan(bbb)], bbb[~np.isnan(bbb)])
plt.title(pol+' ' + field + ' ' + surface + ' ' + canopy + 'R2='+str(r_value)+' RMSE='+str(www))
plt.savefig('/media/wodan/Volume/Arbeit/scatterplot_fertig_'+field+'_'+pol+'_'+file_name+'_'+S.models['surface']+'_'+S.models['canopy'])
pdb.set_trace()



# delete
#----------------------------------------------------------------




high_s = 10*np.log10(pol_field.values.flatten())
high_c = 10*np.log10(S.__dict__['stot'][pol[::-1]])
low_s = 10*np.log10(pol_field.values.flatten())
low_c = 10*np.log10(S.__dict__['stot'][pol[::-1]])
med_s = 10*np.log10(pol_field.values.flatten())
med_c = 10*np.log10(S.__dict__['stot'][pol[::-1]])

aaa = np.append(high_c, low_c)
aaa = np.append(aaa, med_c)
bbb = np.append(high_s, low_s)
bbb = np.append(bbb, med_s)
rm = rmse_prediction(aaa, bbb)

plt.plot(high_s, high_c, 'ks', label='508 High')
plt.plot(low_s, low_c, 'ys', label='508 Low')
plt.plot(med_s, med_c, 'ms', label='508 Medium')
plt.plot(x,x)
plt.xlabel('Backscatter Sentinel-1 [dB]')
plt.ylabel('Backscatter predicted [dB]')
plt.title('Winter Wheat 508, RMSE = 1.17, WCM WCM VH')
plt.legend()
plt.savefig('/media/tweiss/Daten/paper_plot/scatterplot_WCM_WCM_508_vh')




# A_vv = 0.0029
# B_vv = 0.33
# C_vv = [-7.0021433707208622, -6.9312895250907784, -6.5102680381863367, -6.6937605992826743, -6.4430120949713627, -6.0076971888530366, -5.666474504354257, -5.6200494820619369, -5.6494676352946573, -10.100033688395119, -8.9003834449127162, -4.7682144969285822, -4.6781049879584859, -4.2454507330778544, -14.797694536649589, -13.478947387726363, -12.739543885019939, -13.065757732228587, -12.779709894507503, -10.352786544103557, -10.280520800625553, -10.901483934940972, -20.246244459152432, -3.4177109999803292, -3.4027061835701145, -3.4522368335680844, -3.1636592205448153, -2.6195074352905747, -2.6971935752953908, -2.9339919404941948, -4.059268046457758, -4.2960141568979724, -4.5500940632893618, -5.0867665452259283, -5.0698500602677319, -5.7721456271091975, -7.22311510863448, -7.7754890891045951, -6.7785893568690589, -6.3032277546244666, -6.3283266771049282, -6.3884062157800852, -6.3714602721242688, -6.7130767464674559, -4.6582587436532448, -4.6834676736155858, -3.6212693800166251, -2.395594019861313, -1.8295324846832481, -1.9779423468147217, -1.6515616789710934, -0.36079915915934824, 0.010039163680978782, 0.39874983137560466, 0.91329347249675197, 1.6987160760602409, -14.184361799498355, -14.473733017576722, -12.910552531728282, -7.6219662256248535, 3.1365457735983227, 2.9912212384609433, 9.9938626151727341, 10.521504718235484, 13.107744988598153, 17.5773214526007, 8.7424478784475514, 5.3799117334761721]
# D_vv = [-1.0201535399261381, -0.99317394338979437, -0.90633304116716951, -0.94072513784928768, -0.89173568793910363, -0.80786336688157534, -0.75544670018026971, -0.76494878525239096, -0.76067914644358325, 21.034956738407509, 16.448275167290802, -0.61957332700942913, -0.61769417527429038, -0.57352367891223888, 42.140360826181229, 37.045580646946625, 34.534844623512591, 36.173375946129475, 35.573747562752729, 27.691985512512435, 27.47295885766999, 29.789797996818663, 59.969567323536189, -0.56219565623644308, -0.5455222460853707, -0.55402206100894691, -0.48172490562717873, -0.32553249272117457, -0.33684259866002492, -0.37735333144577377, -0.65056441482222382, -0.66642725940592407, -0.73262998349842456, -0.85498046546635131, -0.82483164658902297, -1.0423196399310231, -1.420158229877065, -1.5645677408997394, -1.3176573014748683, -1.1152122470245696, -1.1143985793948521, -1.1470732664940302, -1.1465893335545905, -1.2055024265100753, -0.8081741430996342, -0.79770794599536643, -0.55001065599609633, -0.26836747465739941, -0.20699802459703379, -0.1884770839859996, -0.033586630964986756, 0.29745818205128116, 0.40000654910153888, 0.43299944999683032, 0.56959784176926698, 0.81996272418179939, 71.103566710274066, 69.075840839555966, 62.886028237081561, 43.533314090334962, 1.1510406981929118, 1.1167131312188456, -25.795135128215747, -25.07005915743321, -35.11432361714008, -55.825209642877404, -14.706288420383506, 1.6201773958005654]
# # surface = 'Oh92'
# # surface = 'Oh04'
# # surface = 'Dubois95'
# surface = 'WaterCloud'
# # canopy = 'turbid_isotropic'
# # canopy = 'turbid_rayleigh'
# canopy = 'water_cloud'
# models = {'surface': surface, 'canopy': canopy}



# A_hv = 0.00080064250078689146
# B_hv = 0.4
# C_hv = [-16.634613310519622, -16.501101769337573, -16.417991624106545, -16.612567776708939, -16.706413635693391, -16.542788562536749, -16.508571565156014, -16.587746206829927, -16.172791601218378, -15.911375094254044, -15.71027401125124, -15.614538717347907, -15.462645388637007, -15.057139947499849, -14.784000075662192, -14.689038450144153, -14.444158016745689, -14.426660376117361, -14.389467831292599, -14.322668476865617, -14.255791190315083, -14.116892656069917, -14.168282640050146, -14.586640020417555, -14.744241608712974, -15.041901580971521, -15.274633100126206, -15.476914136929393, -15.55423461069115, -15.70173657493682, -15.928924135193062, -16.297778785255023, -16.385293604046911, -16.517356500192577, -16.522247242645818, -16.915692587444966, -17.039203913488659, -16.944895956399993, -16.901367245779529, -17.057983931371936, -16.731138209070568, -16.493755710389653, -16.320945913682213, -16.448066985746049, -16.332620399765634, -16.3253677331675, -16.452282389784603, -16.458618271394251, -16.592441491343756, -16.705790793902292, -16.374347043144752, -16.155935182900034, -15.650201151254249, -15.2502952918144, -14.700924171310801, -13.756412483068571, -13.267629488413061, -12.986950203713237, -12.682637908996808, -12.354062659813403, -11.72655724388358, -11.53649783364788, -11.364197950004096, -10.916555519008282, -10.818757006179913, -10.800126480787169, -10.245184850292437, -10.044464832662111]
# D_hv = [-2.9838284371683987, -2.9397726948631853, -2.927306164033411, -2.973815227065808, -2.9980005987139711, -2.9711745929497471, -2.9811124440906576, -3.0479475567001533, -2.9873859955569704, -3.0053729805491169, -3.034307796859975, -3.0986503436850712, -3.1423151618211693, -3.1539084073270525, -3.1987389085502649, -3.3169004915686751, -3.3212141241185194, -3.2911296216411574, -3.2618610395068575, -3.3310924620852047, -3.3643619307538164, -3.4511182610276765, -3.5867397287464233, -3.6629880209157299, -3.6395378604421147, -3.700837859353197, -3.7878127743929455, -3.8355595638433551, -3.8161452571198047, -3.8339672647157141, -3.8697225159074389, -3.9051291427932724, -3.8770638911505029, -3.9005532737122639, -3.8787523190949593, -3.919905896264221, -3.8898149055391147, -3.8733411925847561, -3.8617344951795953, -3.9100920601400579, -3.8171302144101777, -3.9153309260256255, -3.9006602199868672, -3.9503321455414597, -3.9090727840625354, -3.8911383622542055, -3.8807415646955934, -3.8696095625052109, -3.8188601429413795, -3.782234027704241, -3.5334860446920335, -3.4433641702014275, -3.289630301739952, -3.1590621806593093, -2.986102458155655, -2.7012377759978441, -2.5603492039116134, -2.5551931801247867, -2.5248018135811758, -2.4715874620336504, -2.3762178055452434, -2.3235066114283081, -2.2907050519280738, -2.2119028568699899, -2.1719158156952911, -2.1103056134326401, -1.8969028193827557, -1.8080438011792312]





# # oh + water
# A_vv = 0.0029
# B_vv = [0.11520901694978017, 0.13000318815904724, 0.1407397682482841, 0.155796378224481, 0.15289104196611447, 0.14837465941401845, 0.15813682702837462, 0.16589184438894836, 0.17496334528903509, 0.1872586776427482, 0.19449840067231894, 0.19076458546225958, 0.19185595041930015, 0.18775981902434144, 0.18044837989804038, 0.17591320420852841, 0.17638907330999737, 0.18287517259954225, 0.18569034240389459, 0.17901283985140651, 0.18030058081761949, 0.18153229415131944, 0.18986040897579687, 0.21015701539462253, 0.21150195382083281, 0.21847709197603829, 0.21539371344122318, 0.20574861272076808, 0.21135756902171698, 0.22383810912358906, 0.25287981408171306, 0.25956125547393299, 0.26656856013128455, 0.28492542405131172, 0.28282123563538891, 0.29999999999999999, 0.29999999999999999, 0.29999999999999999, 0.29999999999999999, 0.29999999999999999, 0.29999999999999999, 0.29999999999999999, 0.29999999999999999, 0.29999999999999999, 0.29999999999999999, 0.29999999999999999, 0.26302015514753874, 0.23834742160688421, 0.22668918871299845, 0.22636090147937826, 0.21002761158607092, 0.18248684307193355, 0.16873024967158792, 0.16490093621273544, 0.1487351361849063, 0.1302196759428563, 0.11746998189525498, 0.12295013722738055, 0.1167172460670486, 0.1041547108118569, 0.09573519362479041, 0.096078486282794154, 0.080173596340211706, 0.066139876003726741, 0.058661907024993874, 0.057837847311163003, 0.044775939837970416, 0.035233320149668984]
# sss = 0.0095
# surface = 'Oh92'
# # surface = 'Oh04'
# # surface = 'Dubois95'
# # surface = 'WaterCloud'
# # canopy = 'turbid_isotropic'
# # canopy = 'turbid_rayleigh'
# canopy = 'water_cloud'
# models = {'surface': surface, 'canopy': canopy}



# A_hv = [0.0024101631029317394, 0.0022414061456099872, 0.0016062812977589013, 0.0011260997320723955, 0.0013433183104701776, 0.0018490634705823396, 0.0050000000000000001, 0.0050000000000000001, 0.0022454949213466619, 0.001875538053930075, 0.0017031400175682686, 0.0017646147675934221, 0.0018390501968671817, 0.0018657936281426288, 0.0018587181529612566, 0.0018511718258667823, 0.0019814335505075479, 0.001843958357716072, 0.0018336394073696623, 0.0018007994314890401, 0.0017173949276130943, 0.001574006140649639, 0.0013773657489671105, 0.001123301251807828, 0.0010088150576047033, 0.00073191484822288523, 0.00065861734816084294, 0.00062181642874059562, 0.00056980633840813643, 0.00046887901712817389, 0.00048524670784220938, 0.00045687177380688513, 0.00042711491122701667, 0.00036067537708343145, 0.0004180486400399193, 0.00039292728696764843, 0.00037675456246068982, 0.00037656340602947513, 0.00042615885437970018, 0.00044916113813820715, 0.00052978518021414695, 0.00050071901464879525, 0.00054812064742052423, 0.00054523935354741627, 0.00056930037232601151, 0.00054520877785682576, 0.00053596959051485529, 0.00055849227911754001, 0.00053348265602644069, 0.00046584692878623885, 0.00060278667838387057, 0.00069586212888543271, 0.00084639278749773202, 0.00094717751931116981, 0.0011714533401233765, 0.0015562835061136968, 0.0018213520618247908, 0.0019659580266260072, 0.0021939365942321301, 0.002460548134067716, 0.0050000000000000001, 0.0050000000000000001, 0.0033402727155514541, 0.0038343246319741332, 0.0039975562372645866, 0.0041037614783891894, 0.0048559764703324515, 0.0050000000000000001]

# B_hv = [0.24999724351092836, 0.24999443747352235, 0.24996856929446762, 0.24991269024015131, 0.24991239174888807, 0.24994590638636768, 0.059608943937142549, 0.055976647788023207, 0.24986545023583476, 0.24990998400401498, 0.24986839470559524, 0.24987779819630418, 0.24989174033588199, 0.24988550981668026, 0.24987145437965466, 0.2498671542369349, 0.24987634052171109, 0.24983538930462693, 0.24984793199462391, 0.24984878291639656, 0.24982445662082661, 0.24976725420013551, 0.2496895229892038, 0.24959960115718624, 0.24953821858995626, 0.24935256123014388, 0.24929117242942572, 0.24929434092047487, 0.24925348085937019, 0.24913293144433468, 0.24913724476093452, 0.24915874638231261, 0.24913367241821605, 0.24903209758395756, 0.24909144806007219, 0.24912588224645518, 0.24910988535252923, 0.24905710587908378, 0.24911533928500673, 0.24917670596620162, 0.24921711593492932, 0.2491469527777502, 0.24920056222738673, 0.24924539555893133, 0.24924709503572401, 0.24917352069368939, 0.24917896605060741, 0.2492352788536519, 0.24920809114844503, 0.24911269250481943, 0.24922345655423819, 0.24931725107638142, 0.24938996646893058, 0.24941441523072602, 0.2495482524088049, 0.24972929218064399, 0.24981158472750334, 0.24984065474208783, 0.24990394153692325, 0.24994931476793822, 0.064326035124241954, 0.072568323821922257, 0.25001491938371884, 0.2501154601930361, 0.2501718770209746, 0.25022756743348923, 0.25069541403761014, 0.35403686834193177]
# s = 0.015




# vv_water = ([-11.29164033, -11.53097481, -10.84022669, -11.87687424,-12.69253232, -12.88291614, -12.22720586, -13.57001433,-14.86199215, -14.87976315, -14.37997695, -15.32978245,-15.20320345, -13.82195555, -14.80531454, -16.04441349,-16.12949497, -14.78256043, -12.64949291, -14.71404135,-15.27605283, -14.10242778, -15.33309722, -16.84346637,-16.5883354 , -15.39481078, -15.93414555, -17.55716044,-17.38124334, -16.30997367, -17.40246767, -18.63149867,-18.42128205, -17.57775692, -17.87494215, -18.9862362 ,-18.94641333, -18.24183742, -18.23065031, -18.92441309,-18.68899449, -17.92580259, -18.11655651, -18.91552272,-18.49185026, -17.64517271, -17.63048388, -18.36609018,-17.89124898, -16.73752631, -16.80911965, -17.57770979,-16.9222695 , -15.14490802, -15.05212064, -16.18656282,-15.96454663, -14.79679952, -14.59644983, -15.05206424,-14.15827177, -12.4317577 , -12.4574119 , -13.60001113,-12.15523458,  -9.83541064, -10.44907831, -12.30497139])


# vv_ssrt = ([-11.73296279, -11.71072028, -10.4099441 , -11.39410103,-13.27592687, -12.99100915, -11.66308522, -12.6834874 ,-14.93613846, -14.88813021, -13.91933916, -15.60611928,-14.86699702, -12.98657047, -13.81063416, -15.67688184,-15.40867325, -13.71464474, -13.76681039, -15.73807507,-15.60025419, -13.84018541, -14.64963066, -17.21830867,-16.82277193, -14.87763574, -15.32866819, -17.85436891,-18.43806798, -17.28605477, -18.89273707, -19.72885073,-19.48595161, -18.14095916, -18.11952393, -19.76923186,-19.71578127, -18.44903697, -18.84344013, -19.96170649,-19.72402981, -18.35878974, -18.53680926, -19.78576501,-19.26198743, -17.40947299, -17.63139079, -19.02971728,-18.60478653, -16.92646248, -16.88307389, -17.98647951,-17.18262845, -15.22711546, -15.073621  , -16.2905793 ,-15.3300122 , -13.69907975, -13.80080936, -14.97275706,-13.67806023, -12.12335029, -12.00410215, -13.31086964,-12.66511182, -11.01034934, -10.80061669, -11.84466433])

# vv_oh = ([-11.85042477, -11.77181586, -10.39360009, -11.21322874,-13.09654463, -12.91849099, -11.579098  , -12.76400063,-15.05283259, -14.9538469 , -13.98244843, -15.69435184,-14.90724022, -12.98956509, -13.87368895, -15.78527483,-15.47075364, -13.72439359, -13.72821759, -15.75256047,-15.56752216, -13.74222141, -14.62093141, -17.27687594,-16.94524728, -15.2021033 , -15.75573244, -17.73707634,-17.570807  , -16.22022987, -17.31112897, -18.89710532,-18.70059062, -17.71003646, -17.99817583, -19.29614209,-19.18259512, -18.57760811, -18.60497818, -19.17100586,-18.94002982, -18.19157125, -18.38259906, -19.1927953 ,-18.79871191, -17.86948892, -17.75673486, -18.60888941,-18.18891001, -16.86886712, -16.85817673, -17.75124838,-17.00199271, -15.18784545, -15.06008892, -16.20218593,-15.27705925, -13.69916554, -13.79799592, -14.94022205,-13.69348997, -12.14686988, -12.02312198, -13.31084253,-12.66951201, -11.0183244 , -10.80230234, -11.84512735])

# fig, ax = plt.subplots(figsize=(20, 10))
# # plt.title('Winter Wheat')
# plt.ylabel('Backscatter [dB]', fontsize=15)
# plt.tick_params(labelsize=12)
# ax.plot(10*np.log10(pol_field), 'ks-', label='Sentinel-1 Pol: VV', linewidth=3)
# ax.plot(date, vv_ssrt, 'C1s-', label='Oh92 + SSRT')
# ax.plot(date, vv_oh, 'rs-', label='Oh92 + Water Cloud (Vegetation)')
# ax.plot(date, vv_water, 'cs-', label='Water Cloud (Surface + Vegetation)')
# ax.legend()
# ax.legend(loc=2, fontsize=12)
# ax.grid(linestyle='-', linewidth=1)
# ax.grid(b=True, which='minor', linestyle='--', linewidth=0.5)

# days = mdates.DayLocator()
# ax.xaxis.set_minor_locator(days)

# months = MonthLocator()
# ax.xaxis.set_major_locator(months)
# ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b %Y"))
# ax.set_ylim([-22,-8])
# ax.set_xlim(['2017-03-20', '2017-07-15'])
# plt.savefig('/media/tweiss/Daten/plots/plot_vv_combination')



# vh_ssrt = ([-19.67840577, -19.42277156, -18.0874326 , -18.784136  ,-20.49293683, -20.13213538, -18.65765464, -19.30629456,-20.62748915, -20.4997889 , -19.40989391, -20.63636707,-19.62986344, -18.04439177, -19.00206007, -20.65364374,-20.33373564, -18.75712058, -18.21032988, -20.0875454 ,-20.05453204, -18.64203735, -19.24113199, -21.18278517,-20.9481051 , -19.446467  , -19.92907533, -22.14824424,-22.40505421, -21.28537842, -21.97744043, -23.21053266,-22.97402324, -21.78444531, -21.92873777, -23.34982046,-23.13909269, -21.94681763, -22.29488843, -23.60842518,-23.21936016, -21.98751723, -22.18966466, -23.23086875,-22.94321577, -21.33926427, -22.05124494, -23.32148346,-23.13168888, -22.01649487, -21.81360599, -22.99877387,-22.42084838, -20.67813528, -20.47526623, -21.34051821,-20.59703076, -18.94740356, -18.91883699, -19.94917444,-18.53092039, -17.1321589 , -17.45929381, -18.81202117,-18.55865424, -17.02931188, -16.85905051, -18.16399906])
# vh_oh = ([-19.86162116, -19.57112375, -18.13734837, -18.69251071,-20.43000572, -19.99688552, -18.53495906, -19.14043324,-20.49576836, -20.34492092, -19.31183183, -20.53365225,-19.62780736, -18.13564594, -19.00995295, -20.52169587,-20.12733427, -18.73005616, -18.34945469, -20.09651229,-20.03063923, -18.67878957, -19.2757966 , -21.29429226,-21.14923855, -19.68089443, -20.24810453, -22.58056418,-22.45110578, -20.81696221, -21.12618394, -23.36930194,-23.21050992, -21.51809928, -21.84881146, -23.92022565,-23.55360833, -21.66719502, -22.24945342, -24.07931536,-23.36632875, -21.63576534, -22.0153772 , -23.39512932,-23.03313252, -20.89670787, -21.85427929, -23.53243939,-23.24947955, -21.68589889, -21.68717774, -23.08263265,-22.26612662, -20.51709422, -20.44459   , -21.01149744,-20.27311166, -18.9290326 , -18.88919762, -19.53818696,-18.48911009, -17.52960988, -17.56977822, -18.11691683,-17.80145748, -16.87040664, -16.54868817, -17.14364399])
# vh_water = ([-18.70551503, -18.65896705, -18.45558832, -18.96082513,-19.4584029 , -19.4771854 , -19.27048062, -19.79958409,-19.92824465, -19.67185765, -19.40449928, -19.99075842,-19.94072544, -19.13119075, -19.0624501 , -19.64101252,-19.29856842, -18.78607457, -19.26241813, -19.85627206,-19.59426813, -18.93042516, -19.28771661, -20.3407397 ,-20.29032592, -19.96788033, -20.50179834, -21.54404228,-21.40531241, -20.85525185, -21.34432156, -22.2828478 ,-22.07766317, -21.44904713, -21.68349998, -22.64739519,-22.48301502, -21.73052239, -21.93403278, -22.77857646,-22.40818623, -21.60901865, -21.71591686, -22.6346061 ,-22.29920863, -21.65404604, -21.77309582, -22.53463294,-22.31702525, -21.6228732 , -21.5903788 , -22.2278543 ,-21.66509546, -20.64525213, -20.25865423, -20.55611041,-19.89034581, -18.86794642, -18.79026734, -19.39460539,-18.69921614, -17.66310926, -17.64455483, -18.06052091,-17.64990205, -16.80258795, -16.43289511, -17.08132995])

# fig, ax = plt.subplots(figsize=(20, 10))
# # plt.title('Winter Wheat')
# plt.ylabel('Backscatter [dB]', fontsize=15)
# plt.tick_params(labelsize=12)
# ax.plot(10*np.log10(pol_field), 'ks-', label='Sentinel-1 Pol: VH', linewidth=3)
# ax.plot(date, vh_ssrt, 'C1s-', label='Oh92 + SSRT')
# ax.plot(date, vh_oh, 'rs-', label='Oh92 + Water Cloud (Vegetation)')
# ax.plot(date, vh_water, 'cs-', label='Water Cloud (Surface + Vegetation)')
# ax.legend()
# ax.legend(loc=2, fontsize=12)
# ax.grid(linestyle='-', linewidth=1)
# ax.grid(b=True, which='minor', linestyle='--', linewidth=0.5)

# days = mdates.DayLocator()
# ax.xaxis.set_minor_locator(days)

# months = MonthLocator()
# ax.xaxis.set_major_locator(months)
# ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b %Y"))
# ax.set_ylim([-25,-15])
# ax.set_xlim(['2017-03-20', '2017-07-15'])
# plt.savefig('/media/tweiss/Daten/plots/plot_vh_combination')


