from astropy import units as u

import numpy as np

import uttr


@uttr.s(repr=False)
class SphericalVoids:
    """
    Represents a collection of spherical voids with associated properties.

    Attributes
    ----------
    rad : numpy.ndarray
        Radius of each void in megaparsecs (Mpc).
    x_void, y_void, z_void : numpy.ndarray
        Center coordinates of each void in Mpc.
    vel_x_void, vel_y_void, vel_z_void : numpy.ndarray
        Velocity components of each void in Mpc/h.
    delta : numpy.ndarray
        Density contrast of each void.
    dtype : numpy.ndarray
        Data type of the input arrays.
    poisson : numpy.ndarray
        Poisson noise for each void.
    nran : numpy.ndarray
        Number of random points used for each void.

    Methods
    -------
    __len__()
        Returns the number of voids in the collection.

    __repr__()
        Returns a string representation of the object."""

    rad = uttr.ib(converter=np.array, unit=u.Mpc)
    x_void = uttr.ib(converter=np.array, unit=u.Mpc)
    y_void = uttr.ib(converter=np.array, unit=u.Mpc)
    z_void = uttr.ib(converter=np.array, unit=u.Mpc)
    vel_x_void = uttr.ib(converter=np.array, unit=u.Mpc / u.h)
    vel_y_void = uttr.ib(converter=np.array, unit=u.Mpc / u.h)
    vel_z_void = uttr.ib(converter=np.array, unit=u.Mpc / u.h)
    delta = uttr.ib(converter=np.array)
    dtype = uttr.ib(converter=np.array)
    poisson = uttr.ib(converter=np.array)
    # dist4 = uttr.ib(converter=np.array)
    nran = uttr.ib(converter=np.array)
    _void_len = uttr.ib(init=False)

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        lengths = set(())
        for e in (
            self.rad,
            self.x_void,
            self.y_void,
            self.z_void,
            self.vel_x_void,
            self.vel_y_void,
            self.vel_z_void,
            self.delta,
            self.dtype,
            self.poisson,
            # self.dist4,
            self.nran,
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
