"""
This module implements the Oh et al. (1992)
empirical surface backscattering model as documented
in Ulaby et al (2014), chapter 10.5

References
----------
Oh et al. (1992): An empirical model and an inversion technique for radar scattering from bare soil surfaces. IEEE TGRS 30(2). 370-381.
"""

import numpy as np
import matplotlib.pyplot as plt

from .. core import Fresnel0
from .. core import Reflectivity

class Oh92(object):
    def __init__(self, eps, ks, theta):
        """
        Parameters
        ----------
        eps : complex
            relative dielectric permitivity
        ks : float
            product of wavenumber and rms height
            be aware that both need to have the same units
        theta : float, ndarray
            incidence angle [rad]
        """

        self.eps = eps
        self.ks = ks
        self.theta = theta

        # calculate p and q
        self.G0 = Fresnel0(self.eps)  # nadir fresnel reflectivity
        self.G = Reflectivity(self.eps, self.theta)
        self._calc_p()
        self._calc_q()

        # calculate backascatter
        self.vv = self._calc_vv()
        self.hh = self.p * self.vv
        self.hv = self.q * self.vv

    def _calc_p(self):
        a = 1./(3.*self.G0.x)
        self.p = (1.- (2.*self.theta/np.pi)**a * np.exp(-self.ks)     )**2.

    def _calc_q(self):
        self.q = 0.23*(self.G0.x)**0.5 * (1.-np.exp(self.ks))

    def _calc_vv(self):

        a = 0.7*(1.-np.exp(-0.65*self.ks**1.8))
        b = np.cos(self.theta)**3. * (self.G.v+self.G.h) / np.sqrt(self.p) 
        return a*b

    def plot(self):
        f = plt.figure()
        ax = f.add_subplot(111)
        t = np.rad2deg(self.theta)
        ax.plot(t, self.hh, color='blue', label='hh')
        ax.plot(t, self.vv, color='red', label='vv')
        ax.plot(t, self.hv, color='green', label='hv')
        ax.grid()
        ax.legend()





