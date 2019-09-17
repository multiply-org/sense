
import pdb
import gp_emulator

import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import sys
import pickle

import os
import pandas as pd

import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator

from sense.canopy import OneLayer
from sense.soil import Soil
from sense import model
import scipy.stats
from scipy.optimize import minimize, root

import datetime



clay = 0.08
sand = 0.12
bulk = 1.5
f = 5.405



"""-----------------------------------
Save and load settings
-----------------------------------"""


def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


"""-----------------------------------
Plot settings
-----------------------------------"""


def plot_config():
    """Update the MPL configuration"""
    config_json = '''{
            "lines.linewidth": 2.0,
            "axes.edgecolor": "#bcbcbc",
            "patch.linewidth": 0.5,
            "legend.fancybox": true,
            "axes.prop_cycle": "cycler('color', ['#FC8D62','#66C2A5', '#8DA0CB','#E78AC3','#A6D854','#FFD92F','#E5C494','#B3B3B3'])",
            "axes.facecolor": "w",
            "axes.labelsize": "large",
            "axes.grid": false,
            "patch.edgecolor": "#eeeeee",
            "axes.titlesize": "x-large",
            "svg.fonttype": "path",
            "xtick.direction" : "out",
            "ytick.direction" : "out",
            "xtick.color": "#262626",
            "ytick.color": "#262626",
            "axes.edgecolor": "#262626",
            "axes.labelcolor": "#262626",
            "axes.labelsize": 12,
            "font.size": 12,
            "legend.fontsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12

    }
    '''
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['xtick.minor.size'] = 10
    plt.rcParams['xtick.minor.width'] = 0.5
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['ytick.minor.size'] = 10
    plt.rcParams['ytick.minor.width'] = 0.5
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica']

    s = json.loads(config_json)
    plt.rcParams.update(s)
    plt.rcParams["axes.formatter.limits"] = [-5, 5]


def pretty_axes(ax):
    """This function takes an axis object ``ax``, and makes it purrty.
    Namely, it removes top and left axis & puts the ticks at the
    bottom and the left"""

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    loc = plt.MaxNLocator(6)
    ax.yaxis.set_major_locator(loc)
    ax.xaxis.set_major_locator(loc)

    ax.tick_params(axis="both", which="both", bottom="on", top="off",
                   labelbottom="on", left="on", right="off", labelleft="on")


plot_config()


"""-----------------------------------
RT-Model
-----------------------------------"""


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

# def wrapper(x):

#     sm = x[0][0]
#     LAI = x[0][1]
#     CAB = x[0][2]
#     theta = np.deg2rad(x[0][3])

#     s = 0.02706007
#     m_vv = 0.68963347
#     omega_vv = 0.02892507
#     m_vh = 0.32105396
#     omega_vh = 0.00370084

#     coef_vv = m_vv * CAB
#     coef_vh = m_vh * CAB
#     SAR = [f, theta]

#     surface_vv    = [clay, sand, bulk, sm, s]
#     vegetation_vv = [1, LAI, coef_vv, omega_vv]


#     surface_vh    = [clay, sand, bulk, sm, s]
#     vegetation_vh = [1, LAI, coef_vh, omega_vh]

#     vv = run_SSRT(surface_vv, vegetation_vv, SAR)[0]

#     vh = run_SSRT(surface_vh, vegetation_vh, SAR)[1]

#     return 10*np.log10(vv)

def wrapper_2(sm,LAI,CAB,theta):

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

    # # r=field
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

    # return 10*np.log10(vv), 10*np.log10(vh), np.rad2deg(theta)
    # return 10*np.log10(vh), LAI, CAB, np.rad2deg(theta)
    return 10*np.log10(vv), 10*np.log10(vh), LAI, CAB, np.rad2deg(theta)
    # return 10*np.log10(vv), LAI, CAB, np.rad2deg(theta)


# pdb.set_trace()



"""-----------------------------------
Emulator
-----------------------------------"""

parameters = ['sm', 'LAI', 'CAB', 'theta']
min_vals = [0.01, 0, 0, 30]
max_vals = [0.55, 7, 1, 45]

n_train = 300
n_validate = 1000

# retval = []
# for iband in range(1):
#     x = gp_emulator.create_emulator_validation(
#         wrapper, parameters, min_vals, max_vals, n_train, n_validate, do_gradient=True, n_tries=15, n_procs=2)
#     retval.append(x)

# gp, validate, validate_output, validate_gradient, emulated_validation, emulated_gradient = x

# ymin = -25
# fig1, axs1 = plt.subplots(nrows=1, ncols=1, figsize=(14, 7))

# slope, intercept, r_value, p_value, std_err = linregress(validate_output, emulated_validation.squeeze())
# axs1.plot(validate_output, emulated_validation, 'o',
#           mec="#FC8D62", mfc="none", rasterized=True)
# axs1.plot([ymin, 0], [ymin, 0], 'k--', lw=0.5)

# p = np.polyfit(validate_output, emulated_validation, 1)
# mae = np.abs(validate_output - emulated_validation.squeeze()).max()
# print("%d & %6.3f & %6.3f & %6.3f & %6.3e & %6.3e\\\\" % (1, slope, intercept, r_value, std_err, mae),)
# x = np.linspace(ymin, 0, 5)
# axs1.plot(x, np.polyval(p, x), '-', lw=0.4)
# axs1.set_ylim(ymin, 0)
# axs1.set_xlim(ymin, 0)
# pretty_axes(axs1)
# axs1.set_title("SSRT_S2")
# plt.tight_layout()

# fig1.savefig("emulator_ssrt_S2.pdf", dpi=600,
#              rasterize=True, bbox_inches="tight")




# gp.save_emulator('emulator_ssrt_s2')



# Create the training samples
training_samples, distributions = gp_emulator.create_training_set(parameters, min_vals, max_vals, n_train=n_train)
# Create the validation samples
validation_samples = gp_emulator.create_validation_set(distributions, n_validate=n_validate)


# Generate the reflectance training set by running the RT model
# for each entry in the training set
training_s1 = np.zeros((n_train, 5))
# training_s1 = np.zeros((n_train, 4))
# training_s1 = np.zeros((n_train, 3))
for i, p in enumerate(training_samples):
    training_s1[i, :] = wrapper_2(p[0], p[1], p[2], p[3])


# Generate the reflectance validation set by running the RT model
# for each entry in the validation set
validation_s1 = np.zeros((n_validate, 5))
# validation_s1 = np.zeros((n_validate, 4))
# validation_s1 = np.zeros((n_validate, 3))
for i, p in enumerate(validation_samples):
    validation_s1[i, :] = wrapper_2(p[0], p[1], p[2], p[3])

# Define and train the emulator from reflectance to LAI
gp = gp_emulator.GaussianProcess(inputs=training_s1, targets=training_samples[:, 0])
gp.learn_hyperparameters(n_tries=15, verbose=False)

# Predict the LAI from the reflectance
ypred, _, _ = gp.predict(validation_s1)

# Plot
fig = plt.figure(figsize=(7,7))
plt.plot(validation_samples[:, 0], ypred, 'o', mfc="none")
plt.plot([min_vals[0], max_vals[0]], [min_vals[0], max_vals[0]],
        '--', lw=3)
x = np.linspace(min_vals[0], max_vals[0], 100)

regress = scipy.stats.linregress(validation_samples[:, 0], ypred)
plt.plot(x, regress.slope*x + regress.intercept, '-')
# plt.xlabel(r"Validation LAI $[m^{2}m^{-2}]$")
# plt.ylabel(r"Retrieved LAI $[m^{2}m^{-2}]$")
plt.title("Slope=%8.4f, "%(regress.slope) +
          "Intercept=%8.4f, "%(regress.intercept) +
          "$R^2$=%8.3f" % (regress.rvalue**2))

gp.save_emulator('/media/tweiss/Daten/buffer_2/emulator_ssrt_s2_rfieldbuffer50')
