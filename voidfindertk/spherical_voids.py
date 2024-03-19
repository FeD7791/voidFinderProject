from astropy import units as u

import numpy as np

import uttr


@uttr.s(repr = False) 
class SphericalVoids:
    rad = uttr.ib(converter=np.array,unit=u.Mpc)
    x_void = uttr.ib(converter=np.array,unit=u.Mpc)
    y_void = uttr.ib(converter=np.array,unit=u.Mpc)
    z_void = uttr.ib(converter=np.array,unit=u.Mpc)
    vel_x_void = uttr.ib(converter=np.array, unit=u.Mpc / u.h)
    vel_y_void = uttr.ib(converter=np.array, unit=u.Mpc / u.h)
    vel_z_void = uttr.ib(converter=np.array, unit=u.Mpc / u.h)
    delta = uttr.ib(converter=np.array)
    dtype = uttr.ib(converter=np.array)
    poisson = uttr.ib(converter=np.array)
    #dist4 = uttr.ib(converter=np.array)
    nran = uttr.ib(converter=np.array)
    _void_len = uttr.ib(init=False)

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        lengths = set(())
        for e in (
            self.rad, 
            self.x_void, self.y_void, 
            self.z_void, self.vel_x_void, self.vel_y_void, self.vel_z_void,
            self.delta,
            self.dtype,
            self.poisson,
            #self.dist4,
            self.nran
            ):
            lengths.add(len(e))

        if len(lengths) != 1:
            raise ValueError("Arrays should be of the same size")

        super().__setattr__("_void_len", lengths.pop())

    def __len__(self):
        """Length method.

        Returns
        -------
            int
                the number of elements in SphericalVoids
        """
        return self._void_len

    def __repr__(self):
        """Representation method.

        Returns
        -------
            str
                Name plus number of points in the box
        """
        cls_name = type(self).__name__
        length = len(self)
        return f"<{cls_name} size={length}>"
