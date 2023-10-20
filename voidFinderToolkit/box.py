import numpy as np
from astropy import units as u
import uttr
from ..Tools.box_tools import slicer


@uttr.s
class Box:
    """
    Class used to describe a set of points (x,y,z) alongside with its
    velocities (vx,vy,vz)

    Attributes
    ----------
    x : numpy.ndarray
    y : numpy.ndarray
    z : numpy.ndarray
        (x,y,z) array of position elements
    vx : numpy.ndarray
    vy : numpy.ndarray
    vz : numpy.ndarray
        (vx,vy,vz) array of velocity elements
    m : numpy.ndarray

    Methods
    -------

    slice(self,n,parameter)
        Calls the function Slicer to slice the set in n parts
    """

    x = uttr.ib(converter=np.array, unit=u.Mpc)
    y = uttr.ib(converter=np.array, unit=u.Mpc)
    z = uttr.ib(converter=np.array, unit=u.Mpc)
    vx = uttr.ib(converter=np.array, unit=u.Mpc / u.h)
    vy = uttr.ib(converter=np.array, unit=u.Mpc / u.h)
    vz = uttr.ib(converter=np.array, unit=u.Mpc / u.h)
    m = uttr.ib(converter=np.array, unit=u.M_sun)

    def __attrs_post_init__(self):
        try:
            length_set = {
                len(self.x),
                len(self.y),
                len(self.z),
                len(self.vx),
                len(self.vy),
                len(self.vz),
                len(self.m),
            }
            if len(length_set) != 1:
                raise ValueError('Arrays should have the same size')
        except ValueError as error:
            print(error)
            print('Arrays sizes were:',
            len(self.x),
            len(self.y),
            len(self.z),
            len(self.vx),
            len(self.vy),
            len(self.vz),
            len(self.m))

    def slice(self, n, parameter):
        """
        Calls the function Slicer to divide the dataset in n parts by
        classifying the points based on the distance to the origin
        (0,0,0)

        Parameters
        ----------

        n : int
            Number of divitions of each axis
        parameter : str {'position','velocity'}
            'position' divides the dataset based on the (x,y,z)
            coordinates of the dot cloud
            'velocity' divides the dataset based on the (vx,vy,vz)
            coordinates of the dot cloud
        """
        if parameter == "position":
            columns = ["x", "z", "y", "distance"]
            slices = slicer(
                self.x / u.Mpc, self.y / u.Mpc, 
                self.z / u.Mpc, n, columns
                )
            slices = map(np.array, slices)
            slices = map(lambda a: a * u.Mpc, slices)
        elif parameter == "velocity":
            columns = ["vx", "vz", "vy", "v_modulus"]
            slices = slicer(
                self.vx * u.h / u.Mpc,
                self.vy * u.h / u.Mpc,
                self.vz * u.h / u.Mpc,
                n,
                columns,
            )
            slices = map(np.array, slices)
            slices = map(lambda a: a * u.Mpc / u.h, slices)
        return list(slices)
