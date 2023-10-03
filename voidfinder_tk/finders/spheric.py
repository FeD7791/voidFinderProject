# finders/spheric.py
from typing import NamedTuple

class SphericalVoid(NamedTuple):
    radius : float
    x : float
    y : float
    z : float
    vx : float
    vy : float
    vz : float
    delta : float
    d_type : float
    
    def type_of_v(self):
        delta = self.d_type
        _type_void : str
        try:
            if (delta>0):
                _type_void = 'S'
            elif (delta<0):
                _type_void = 'R'
            else:
                return 'No type'
        except:
            print('Invalid delta')
        else:
            return _type_void
