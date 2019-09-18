
import pandas as pd
import os
import pdb
import matplotlib.pyplot as plt
import gdal
import kafka
import numpy as np
import fnmatch
import subprocess
import datetime
from kafka.input_output import Sentinel1_Observations
import gp_emulator
import matplotlib
import glob

from multiply_prior_engine.interpolate_grids import interpolate_stack
from multiply_prior_engine.interpolate_grids import get_files

from sense.canopy import OneLayer
from sense.soil import Soil
from sense import model

import scipy.stats

def _error(actual: np.ndarray, predicted: np.ndarray):
    """ Simple error """
    return actual - predicted

def me(actual: np.ndarray, predicted: np.ndarray):
    """ Mean Error """
    return np.nanmean(_error(actual, predicted))

def mse(actual: np.ndarray, predicted: np.ndarray):
    """ Mean Squared Error """
    return np.nanmean(np.square(_error(actual, predicted)))

def rmse(actual: np.ndarray, predicted: np.ndarray):
    """ Root Mean Squared Error """
    return np.sqrt(mse(actual, predicted))

def ubrmse(actual: np.ndarray, predicted: np.ndarray):
    """ unbiased Root Mean Squared Error """
    return np.sqrt(np.nansum((actual-(predicted-np.nanmean(actual)))**2)/len(actual))

def linregress(predictions, targets):
    """ Calculate a linear least-squares regression for two sets of measurements """
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(predictions, targets)
    return slope, intercept, r_value, p_value, std_err

def filter_relativorbit(data, field, orbit1, orbit2=None, orbit3=None, orbit4=None):
    """ data filter for relativ orbits """
    output = data[[(check == orbit1 or check == orbit2 or check == orbit3 or check == orbit4) for check in data[(field,'relativeorbit')]]]
    return output

def list_files(inputfolder, expression=None):
    """Create list containing all files in inputfolder (without subfolder) that contain expression in the name"""
    if expression is None:
        expression = '*'
    filelist = glob.glob(os.path.join(inputfolder,expression))
    print("Number of found files:", len(filelist))
    return filelist

def run_SSRT(surface, vegetation, SAR):
    """
    implementation of Oh model and Single Scattering RT-model
    Return Polarization VV and VH in linear scale
    """

    clay, sand, bulk, sm, s = surface
    d, LAI, coef, omega = vegetation
    f, theta = SAR

    ke = coef*np.sqrt(LAI)

    soil = Soil(mv=sm, s=s, clay=clay, sand=sand, f=f, bulk=bulk)

    can = OneLayer(canopy='turbid_isotropic', ke_h=ke, ke_v=ke, d=d, ks_h = omega*ke, ks_v = omega*ke)

    S = model.RTModel(surface=soil, canopy=can, models= {'surface': 'Oh92', 'canopy': 'turbid_isotropic'}, theta=theta, freq=f)
    S.sigma0()

    vv = S.__dict__['stot']['vv']
    vh = S.__dict__['stot']['hv']

    return vv, vh

def wrapper_2(sm,LAI,CAB,theta):
    """
    wrapper around run_SSRT
    implementation of empirical parameters

    Return:
    - VV (dB), ndarray
    - VH (dB), ndarray
    - LAI, ndarray
    - CAB, ndarray
    - theta (degree), ndarray
    """

    theta = np.deg2rad(theta)

    # s = 0.02706007
    # m_vv = 0.68963347
    # omega_vv = 0.02892507
    # m_vh = 0.32105396
    # omega_vh = 0.00370084

    # #r=pixel (1x1)
    # s = 0.03
    # m_vv = 1.17135053
    # omega_vv = 0.04116796
    # m_vh = 0.56848851
    # omega_vh = 0.00857545

    # # r=15m (4x4)
    # s = 0.03
    # m_vv = 1.20086761
    # omega_vv = 0.04472383
    # m_vh = 0.63712893
    # omega_vh = 0.01117671

    # # r=50m (?x?)
    # s = 0.03
    # m_vv = 1.29354856
    # omega_vv = 0.05000378
    # m_vh = 0.57818271
    # omega_vh = 0.0094646

    # # r=field
    # s = 0.03
    # m_vv = 1.40816754
    # omega_vv = 0.07
    # m_vh = 0.84442888
    # omega_vh = 0.02004628

    # # r=field_CAB
    # s = 0.01758208
    # m_vv = 2.84339935
    # omega_vv = 0.07
    # m_vh = 1.19354843
    # omega_vh = 0.01752458

    # # r=fieldbuffer50
    # s = 0.03
    # m_vv = 1.35998835
    # omega_vv = 0.05649291
    # m_vh = 0.61930753
    # omega_vh = 0.01189246

    # r=fieldbuffer50 descending
    s = 0.03
    m_vv = 1.29791155
    omega_vv = 0.05277558
    m_vh = 0.57572605
    omega_vh = 0.01057413

    coef_vv = m_vv * CAB
    coef_vh = m_vh * CAB
    SAR = [f, theta]

    surface_vv    = [clay, sand, bulk, sm, s]
    vegetation_vv = [1, LAI, coef_vv, omega_vv]


    surface_vh    = [clay, sand, bulk, sm, s]
    vegetation_vh = [1, LAI, coef_vh, omega_vh]

    vv = run_SSRT(surface_vv, vegetation_vv, SAR)[0]

    vh = run_SSRT(surface_vh, vegetation_vh, SAR)[1]

    return 10*np.log10(vv), 10*np.log10(vh), LAI, CAB, np.rad2deg(theta)

# Needed parameters for RT-model implementation
clay = 0.08
sand = 0.12
bulk = 1.5
f = 5.405

### Gathering of in-situ measurements (soil moisture, LAI, height, VWC) within Panda dataframe
#------------------------------------------------------------------------------

# LAI path
path = 'data_kafka/insitu/vegetation'

# SM path
path_SM = 'data_kafka/insitu/soil_moisture/'

# field names
# fields = ['301', '508', '542', '319', '515']
fields = ['508']

# ESU names
esus = ['high', 'low', 'med']
esus = ['high']

processed_sentinel_data = 'multi'

df_new = pd.DataFrame()
df_SM = pd.DataFrame(columns=pd.MultiIndex(levels=[[],[]], labels=[[],[]]))
df_SM2 = pd.DataFrame()

for field in fields:
    df_SM = pd.DataFrame(columns=pd.MultiIndex(levels=[[],[]], labels=[[],[]]))

    if field == '301':
        name_SM = {'high': 'Alex4_23May17_17Jul17_SM.csv', 'low': 'Alex5_23May17_17Jul17_SM.csv', 'med': 'Alex6_23May17_17Jul17_SM.csv'}
    if field == '508':
        name_SM = {'high': 'Stud1_23May17_20Jul17_SM.csv', 'low': 'Stud2_23May17_20Jul17_SM.csv', 'med': 'blue_box_2017.csv'}
    if field == '542':
        name_SM = {'high': 'Alex1_23May17_20Jul17_SM.csv', 'low': 'Alex2_23May17_20Jul17_SM.csv', 'med': 'Alex3_23May17_19Jul17_SM.csv'}
    if field == '319':
        name_SM = {'high': 'Philip4_07Jun17_04Sep17_SM.csv', 'low': 'Philip5_07Jun17_04Sep17_SM.csv', 'med': 'Philip6_07Jun17_04Sep17_SM.csv'}
    if field == '515':
        name_SM = {'high': 'Philip1_07Jun17_22Sep17_SM.csv', 'low': 'Philip2_07Jun17_22Sep17_SM.csv', 'med': 'Philip3_07Jun17_22Sep17_SM.csv'}

    df = pd.io.parsers.read_csv(os.path.join(path, field + '.csv'), header=[0, 1], sep=';')
    df = df.set_index(pd.to_datetime(df['None']['None']))
    df_reindexed = df.reindex(pd.date_range(start=df.index.min(), end=df.index.max(), freq='10Min')).interpolate(kind='cubic', axis=0)

    for key in name_SM:
        data_SM = pd.read_csv(os.path.join(path_SM, name_SM[key]))
        data_SM['date'] = data_SM['date'] = pd.to_datetime(data_SM['date'])
        data_SM = data_SM.set_index('date')
        del data_SM.index.name

        if name_SM[key] == 'blue_box_2017.csv':
            data_SM.index = data_SM.index + pd.Timedelta('5 min')
            df_SM[field + ' ' + key,'SM 5cm'] = data_SM[['5TM soil moisture addr 1', '5TM soil moisture addr 2', '5TM soil moisture addr 3']].mean(axis=1)/100.
        else:
            df_SM[field + ' ' + key,'SM 5cm'] = data_SM[['Port1_SM', 'Port2_SM']].mean(axis=1)

        df_add = df_reindexed.filter(like='LAI HANTS')
        df_add = df_add.filter(like=key)
        df_new = df_new.append(df_add)
        df_add = df_reindexed.filter(like='Height [cm]')
        df_add = df_add.filter(like=key)
        df_new = df_new.append(df_add)
        df_add = df_reindexed.filter(like='Water content total kg/m2 HANTS')
        df_add = df_add.filter(like=key)
        df_new = df_new.append(df_add)

        df_SM2 = df_SM2.append(df_SM)

    df_add = df_reindexed.filter(like='LAI mean HANTS')
    df_new = df_new.append(df_add)
    df_add = df_reindexed.filter(like='Height [cm] mean')
    df_new = df_new.append(df_add)
    df_add = df_reindexed.filter(like='Water content total kg/m2 mean HANTS')
    df_new = df_new.append(df_add)

df_new = df_new.append(df_SM2)
df_LAI = df_new.groupby(df_new.index).mean()
#-------------------------------------------------------------------------------

### Run retrieval and gather additional Information
#-------------------------------------------------------------------------------

# path where preprocessd S1-data are stored
path_data = '/media/nas_data/Thomas/S1/processed/MNI_2017/step3'

# path where S2 data for LAI and CAB or stored
S2_path = '/media/tweiss/Daten/UCL_2'

# name part of shapefile (needed spatial aggregation for soil moisture retrieval, set as mask within kafka)
radius = 'fieldbuffer50'

# path of shapefile which is set as mask
path_shapefile = '//media/tweiss/Work/GIT/GitHub/PMarzahn/pygeo/sense/data_kafka/mask_shapefile/'

# name part of shapefile (needed spatial aggregation for soil moisture retrieval, set as mask within kafka)
name_shapefile = 'radius_'

# load gp_emulator for soil moisture retrieval
gp = gp_emulator.GaussianProcess(emulator_file='/media/tweiss/Work/GIT/GitHub/PMarzahn/pygeo/sense/data_kafka/emulator_ssrt_s2_r'+radius+'.npz')

df_output = pd.DataFrame(columns=pd.MultiIndex(index=[], levels=[[],[]], labels=[[],[]]))

for field in fields:

    for esu in esus:

        # rasterize shapefile (mask for calculations)
        subprocess.call('gdal_rasterize -at -of GTiff -burn 1 -te 701514 5347515 701984 5347855 -tr 10 10 -ot Byte -co \"COMPRESS=DEFLATE\" '+path_shapefile+name_shapefile+radius+'m.shp '+path_shapefile+'mask.tif', shell=True)

        name_ESU = 'mask.tif'
        sentinel1_observations = kafka.Sentinel1_Observations.S1Observations(path_data, os.path.join(path_shapefile, name_ESU), emulators={'VH': kafka.observation_operators.sar_observation_operator,'VV': kafka.observation_operators.sar_observation_operator})
        g = gdal.Open(os.path.join(path_shapefile, name_ESU))

        state_mask = g.ReadAsArray().astype(np.bool)

        fl_lai = sorted(list_files(S2_path,'*'+field+'*lai.tif'))
        fl_cab = sorted(list_files(S2_path,'*'+field+'*cab.tif'))

        # reproject S2 data (LAI, CAB)
        for xxx in range(len(fl_lai)):
            S2_lai = Sentinel1_Observations.reproject_image(fl_lai[xxx], os.path.join(path_shapefile, name_ESU))
            arr = S2_lai.ReadAsArray()
            [cols, rows] = arr.shape
            driver = gdal.GetDriverByName("GTiff")
            outdata = driver.Create(fl_lai[xxx][:-4]+'2.tif', rows, cols, 1, gdal.GDT_Float32)
            outdata.SetGeoTransform(g.GetGeoTransform())##sets same geotransform as input
            outdata.SetProjection(g.GetProjection())##sets same projection as input
            outdata.GetRasterBand(1).WriteArray(arr)

            outdata.FlushCache() ##saves to disk!!
            outdata = None
            band=None
            ds=None


            S2_cab = Sentinel1_Observations.reproject_image(fl_cab[xxx], os.path.join(path_shapefile, name_ESU))
            arr = S2_cab.ReadAsArray() / 100.
            [cols, rows] = arr.shape
            driver = gdal.GetDriverByName("GTiff")
            outdata = driver.Create(fl_cab[xxx][:-4]+'2.tif', rows, cols, 1, gdal.GDT_Float32)
            outdata.SetGeoTransform(g.GetGeoTransform())##sets same geotransform as input
            outdata.SetProjection(g.GetProjection())##sets same projection as input
            outdata.GetRasterBand(1).WriteArray(arr)

            outdata.FlushCache() ##saves to disk!!
            outdata = None
            band=None
            ds=None

        # daily interpolation of S2 data for LAI and CAB
        band = 0           # only specify if multiple bands are in tiff however, will be checked anyways)
        fillvalue = -999.  # only specify if multiple bands are in tiff (however, will be checked anyways)

        fl_lai = sorted(list_files(S2_path,'*'+field+'*lai2.tif'))
        stack_lai, dates_lai = interpolate_stack(fl_lai, band, fillvalue)
        fl_cab = sorted(list_files(S2_path,'*'+field+'*cab2.tif'))
        stack_cab, dates_cab = interpolate_stack(fl_cab, band, fillvalue)

        SM = []
        LAI = []
        Height = []
        VWC = []
        LAI2 = []
        Height2 = []
        VWC2 = []
        sigma_sentinel_vv = []
        sigma_sentinel_vh = []
        theta = []
        relativeorbit = []
        orbitdirection = []
        satellite = []
        dates = []

        sm_predict = []

        S2_lai_time_series=[]
        S2_cab_time_series=[]


        for i, tim in enumerate(sorted(sentinel1_observations.dates)):
            if tim < datetime.datetime(2017, 3, 25) or tim > datetime.datetime(2017, 7, 20):
                pass
            else:
                if tim >= df_LAI.index.min() and tim <= df_LAI.index.max():
                    try:
                        data = sentinel1_observations.get_band_data(tim, 'vv_' + processed_sentinel_data)
                        data_vh = sentinel1_observations.get_band_data(tim, 'vh_' + processed_sentinel_data)

                        df_add = df_LAI.filter(like=field).filter(like=esu)


                        LAI.append(df_add.filter(like='LAI HANTS').loc[tim.strftime("%Y-%m-%d %H:00:00")].values[0])
                        SM.append(df_add.filter(like='SM 5cm').loc[tim.strftime("%Y-%m-%d %H:00:00")].values[0])
                        Height.append(df_add.filter(like='Height [cm]').loc[tim.strftime("%Y-%m-%d %H:00:00")].values[0])
                        VWC.append(df_add.filter(like='Water content total kg/m2 HANTS').loc[tim.strftime("%Y-%m-%d %H:00:00")].values[0])

                        df_add2 = df_LAI.filter(like=field)
                        LAI2.append(df_add2.filter(like='LAI mean HANTS').loc[tim.strftime("%Y-%m-%d %H:00:00")].values[0])
                        Height2.append(df_add2.filter(like='Height [cm] mean').loc[tim.strftime("%Y-%m-%d %H:00:00")].values[0])
                        VWC2.append(df_add2.filter(like='Water content total kg/m2 mean HANTS').loc[tim.strftime("%Y-%m-%d %H:00:00")].values[0])
                        print(tim)

                        observations = data.observations*1.
                        observations_vh = data_vh.observations*1.

                        observations[~data.mask] = np.nan
                        observations[~state_mask] = np.nan
                        data.metadata['incidence_angle'][~data.mask] = np.nan
                        data.metadata['incidence_angle'][~state_mask] = np.nan

                        observations_vh[~data.mask] = np.nan
                        observations_vh[~state_mask] = np.nan
                        data_vh.metadata['incidence_angle'][~data.mask] = np.nan
                        data_vh.metadata['incidence_angle'][~state_mask] = np.nan


                        sigma_vv = observations
                        sigma_vh = observations_vh


                        sm_end = 0 * sigma_vv
                        lai_end = np.nan * sigma_vv
                        cab_end = np.nan * sigma_vv


                        S2_lai = stack_lai[dates_lai.index(tim.date())]

                        S2_cab = stack_cab[dates_lai.index(tim.date())]


                        prediction = np.column_stack((10*np.log10(np.nanmean(sigma_vv.flatten())),10*np.log10(np.nanmean(sigma_vh.flatten())), np.nanmean(S2_lai.flatten()),np.nanmean(S2_cab.flatten()),np.nanmean(data.metadata['incidence_angle'].flatten())))

                        ypred, _, _ = gp.predict(prediction)


                        # matplotlib.interactive(False)

                        sm_predict.append(np.nanmean(ypred))

                        S2_lai = stack_lai[dates_lai.index(tim.date())]
                        S2_cab = stack_cab[dates_lai.index(tim.date())]

                        S2_lai[~state_mask] = np.nan
                        S2_cab[~state_mask] = np.nan

                        S2_lai_time_series.append(np.nanmean(S2_lai.flatten()))
                        S2_cab_time_series.append(np.nanmean(S2_cab.flatten()*100))


                        sigma_sentinel_vv.append(np.nanmean(sigma_vv))
                        theta.append(np.nanmean(data.metadata['incidence_angle']))
                        relativeorbit.append(data.metadata['relativeorbit'])
                        orbitdirection.append(data.metadata['orbitdirection'])
                        satellite.append(data.metadata['satellite'])

                        data = sentinel1_observations.get_band_data(tim, 'vh_' + processed_sentinel_data)
                        observations = data.observations*1.
                        observations[~data.mask] = np.nan
                        observations[~state_mask] = np.nan
                        sigma_vh = observations
                        sigma_sentinel_vh.append(np.nanmean(sigma_vh))
                        dates.append(tim)
                    except (KeyError, IndexError):
                        pass

        df_output[field + '_' + esu, 'date'] = dates
        df_output[field + '_' + esu, 'sigma_sentinel_vv'] = sigma_sentinel_vv
        df_output[field + '_' + esu, 'sigma_sentinel_vh'] = sigma_sentinel_vh
        df_output[field + '_' + esu, 'theta'] = theta
        df_output[field + '_' + esu, 'relativeorbit'] = relativeorbit
        df_output[field + '_' + esu, 'orbitdirection'] = orbitdirection
        df_output[field + '_' + esu, 'satellite'] = satellite
        df_output[field + '_' + esu, 'LAI'] = LAI
        df_output[field + '_' + esu, 'SM'] = SM
        df_output[field + '_' + esu, 'Height'] = Height
        df_output[field + '_' + esu, 'VWC'] = VWC
        df_output[field + '_' + esu, 'vh/vv'] = np.array(sigma_sentinel_vh) / np.array(sigma_sentinel_vv)
        df_output[field + '_' + esu, 'S2_LAI'] = S2_lai_time_series
        df_output[field + '_' + esu, 'S2_CAB'] = S2_cab_time_series
        df_output[field + '_' + esu, 'SM_retrieved'] = sm_predict
    df_output[field + '_mean', 'date'] = dates
    df_output[field + '_mean', 'LAI'] = LAI2
    df_output[field + '_mean', 'Height'] = Height2
    df_output[field + '_mean', 'VWC'] = VWC2
    SM2 = df_output.filter(like=field).filter(like='SM').mean(axis=1)
    df_output[field + '_mean', 'SM'] = SM2
    vv = df_output.filter(like=field).filter(like='sigma_sentinel_vv').mean(axis=1)
    vh = df_output.filter(like=field).filter(like='sigma_sentinel_vh').mean(axis=1)
    df_output[field + '_mean', 'sigma_sentinel_vv'] = vv
    df_output[field + '_mean', 'sigma_sentinel_vh'] = vh
    df_output[field + '_mean', 'relativeorbit'] = relativeorbit
    df_output[field + '_mean', 'orbitdirection'] = orbitdirection
    df_output[field + '_mean', 'satellite'] = satellite
    theta_mean = df_output.filter(like=field).filter(like='theta').mean(axis=1)
    df_output[field + '_mean', 'theta'] = theta_mean


save_path = '/media/tweiss/Daten/buffer_2'
df_output.to_csv(os.path.join(save_path, processed_sentinel_data + '_'+radius+'.csv'), encoding='utf-8', sep=',', float_format='%.4f')




### plotting


matplotlib.interactive(False)
fieldname = '508_high'

fig, ax = plt.subplots(figsize=(20, 20))
# plt.title('Winter Wheat')
plt.ylabel('Soil Moisture', fontsize=15)
plt.xlabel('Date', fontsize=15)
plt.tick_params(labelsize=12)
ax.set_ylim([0.15,1])

ax.plot(df_output[fieldname]['date'],df_output[fieldname]['SM'], color='blue', label='soil moisture in-situ')
ax.plot(df_output[fieldname]['date'],df_output[fieldname]['SM_retrieved'], color='orange', label='soil moisture retrieved')


# asm_predict = np.array(sm_predict)
# asm_predict[asm_predict<0] = np.nan
# SM[np.int(np.argwhere(np.isnan(asm_predict)))] = np.nan
# sm_predict[np.int(np.argwhere(np.isnan(asm_predict)))] = np.nan
# SM2 = [incom for incom in SM if str(incom) != 'nan']
# sm_predict2 = [incom for incom in sm_predict if str(incom) != 'nan']


slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.array(df_output[fieldname]['SM']), np.array(df_output[fieldname]['SM_retrieved']))

rmse1 = rmse(np.array(df_output[fieldname]['SM']), np.array(df_output[fieldname]['SM_retrieved']))
ubrmse1 = ubrmse(np.array(df_output[fieldname]['SM']), np.array(df_output[fieldname]['SM_retrieved']))
ubrmse2 = ubrmse(np.array(df_output[fieldname]['SM'][0:29]), np.array(df_output[fieldname]['SM_retrieved'][0:29]))
rmse2 = rmse(np.array(df_output[fieldname]['SM'][0:29]), np.array(df_output[fieldname]['SM_retrieved'][0:29]))

ax.set_title('Soil Moisture: RMSE:'+str(rmse1)+', ubRMSE:'+str(ubrmse1)+', R2:'+str(r_value)+',rmse first halb:'+str(rmse2)+' ubrmse first half:'+str(ubrmse2),fontsize=15)
ax1 = ax.twinx()

ax1.plot(df_output[fieldname]['date'],10*np.log10(df_output[fieldname]['sigma_sentinel_vv']), color='black')
ax1.plot(df_output[fieldname]['date'],10*np.log10(np.array(df_output[fieldname]['sigma_sentinel_vh'])), color='black', label = 'S1 data')


vv, vh, LAI, CAB , theta = wrapper_2(np.array(df_output[fieldname]['SM']), np.array(df_output[fieldname]['S2_LAI']), np.array(df_output[fieldname]['S2_CAB'])/100, np.array(df_output[fieldname]['theta']))

ax1.plot(df_output[fieldname]['date'],vv, color='red')
ax1.plot(df_output[fieldname]['date'],vh, color='red', label='model output when using soil moisture as input')
ax1.set_ylim([-35,-10])
ax.grid()
ax2 = ax.twinx()
ax2.plot(df_output[fieldname]['date'],np.array(df_output[fieldname]['S2_LAI']),color='green', label='LAI')
ax2.set_ylim([0,5])
ax3 = ax.twinx()
ax3.plot(df_output[fieldname]['date'],np.array(df_output[fieldname]['S2_CAB'])/100,color='yellow', label='CAB')
ax3.set_ylim([0,1])

ax.legend()
ax1.legend()
ax2.legend()
ax3.legend()

plt.savefig('/media/tweiss/Work/GIT/GitHub/PMarzahn/pygeo/sense/data_kafka/soil_moisture_retrieval_r'+radius+'m.png')
