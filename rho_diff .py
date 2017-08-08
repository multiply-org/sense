#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 15:34:04 2017

@author: bmueller
"""
import numpy as np


eps1 = 1.
eps2 = 3.-0.j
theta1 = np.pi/6

#from ReflTransm

sin_theta2 = np.sqrt(eps1)/np.sqrt(eps2)*np.sin(theta1)
cos_theta2 = np.sqrt(1 - sin_theta2**2)

rhoh = (np.sqrt(eps1)*np.cos(theta1)-np.sqrt(eps2)*cos_theta2) / \
    (np.sqrt(eps1)*np.cos(theta1) + np.sqrt(eps2)*cos_theta2)
rhov = (np.sqrt(eps1)*cos_theta2-np.sqrt(eps2)*np.cos(theta1)) / \
    (np.sqrt(eps1)*cos_theta2 + np.sqrt(eps2)*np.cos(theta1))

print(rhoh)
print(rhov)

#from Prism1

n1 = np.sqrt(eps1)
n2 = np.sqrt(eps2)
costh2 = np.sqrt(1 - (n1*np.sin(theta1)/n2)**2)


rho_v = -(n2*np.cos(theta1) - n1*costh2)/(n2*np.cos(theta1)+n1*costh2)
rho_h = (n1*np.cos(theta1) - n2*costh2)/(n1*np.cos(theta1)+n2*costh2)

print(rho_h)
print(rho_v)
