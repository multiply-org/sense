"""
Major surface scatter class
"""
class SurfaceScatter(object):
    def __init__(self, eps=None, ks=None, theta=None, kl=None, mv=None, C_hh=None, C_vv=None, D_hh=None, D_vv=None, **kwargs):
        self.eps = eps
        self.ks = ks
        self.theta = theta
        self.kl = kl

        self.mv = mv
        self.C_hh = C_hh
        self.C_vv = C_vv
        self.D_hh = D_hh
        self.D_vv = D_vv

        self._check()

    def _check(self):
        pass
        # assert isinstance(self.eps, complex)

