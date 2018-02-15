"""
Major surface scatter class
"""
class SurfaceScatter(object):
    def __init__(self, eps, ks, theta, kl=None, **kwargs):
        self.eps = eps
        self.ks = ks
        self.theta = theta
        self.kl = kl
        self._check()

    def _check(self):
        assert isinstance(self.eps[0], complex)

class SurfaceScatterWaterCloud(object):
    def __init__(self, mv, theta, C_hh=None, C_vv=None, D_hh=None, D_vv=None):
        self.mv = mv
        self.theta = theta
        self.C_hh = C_hh
        self.C_vv = C_vv
        self.D_hh = D_hh
        self.D_vv = D_vv

    def _check(self):
        #needs to be implemented
        pass



