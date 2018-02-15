"""
Specification of canopies
"""

class Canopy(object):
    def __init__(self, **kwargs):
        self.d = kwargs.get('d', None)
        self._check()

    def _check(self):
        assert self.d is not None, 'Vegetation height needs to be given'



class OneLayer(Canopy):
    """
    define a homogeneous one layer canopy
    """
    def __init__(self, **kwargs):
        super(OneLayer, self).__init__(**kwargs)
        self.ke_h = kwargs.get('ke_h', None)
        self.ke_v = kwargs.get('ke_v', None)
        assert self.ke_h is not None
        assert self.ke_v is not None
        self.ks_h = kwargs.get('ks_h', None)
        self.ks_v = kwargs.get('ks_v', None)
        assert self.ks_h is not None
        assert self.ks_v is not None


class WaterCloudLayer(object):
    """
    """
    def __init__(self, **kwargs):
        self.A_hh = kwargs.get('A_hh', None)
        self.B_hh = kwargs.get('B_hh', None)
        assert self.A_hh is not None
        assert self.B_hh is not None
        self.A_vv = kwargs.get('A_vv', None)
        self.B_vv = kwargs.get('B_vv', None)
        assert self.A_vv is not None
        assert self.B_vv is not None
        self.V1 = kwargs.get('V1', None)
        self.V2 = kwargs.get('V2', None)
        assert self.V1 is not None
        assert self.V2 is not None
