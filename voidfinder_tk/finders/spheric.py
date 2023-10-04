from typing import NamedTuple


class SphericalVoid(NamedTuple):

    '''
    Class definition for a Sperical Void Object

    Parameters
    ----------
    radius : float
        Radius of the Spherical Void
    x : scalar(float)
        x position of the center of the Spherical Void
    y : scalar(float)
        y position of the center of the Spherical Void
    z : scalar(float)
        z position of the center of the Spherical Void
    vx : scalar(float)
        medium velocity of the void in the x direction
    vy : scalar(float)
        medium velocity of the void in the y direction
    vz : scalar(float)
        medium velocity of the void in the z direction
    d_type : scalar(float)
        density for radius > radius void
    '''
    radius: float
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    d_type: float

    def type_of_v(self):
        '''
        Function used to calculate the type of the void. Types could
        be: S (High density profile) , R (Low density profile).

        Returns
        -------
        type_void : string
            S : High density profile void
            R : Low density profile void


        '''  
        _delta = self.d_type
        type_void: str
        try:
            if (_delta > 0):
                type_void = 'S'
            elif (_delta < 0):
                type_void = 'R'
            else:
                return 'No type'   
        except (TypeError):
            print('Invalid delta')
        else:
            return type_void
