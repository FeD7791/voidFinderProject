
import uttr

from astropy import units as u

@uttr.s(repr=False, frozen=True)
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

    Methods
    -------
    build_df(self)
        Transforms the set into a pandas dataframe

    _check_length(self)
        Verifies that the lenght of the inputs are the same

    slice(self,n,parameter)
        Calls the function Slicer to slice the set in n parts
    """

    x = uttr.ib(unit=u.Mpc)
    y = uttr.ib(unit=u.Mpc)
    z = uttr.ib(unit=u.Mpc)
    vx = uttr.ib(unit=u.Mpc/u.h)
    vy = uttr.ib(unit=u.Mpc/u.h)
    vz = uttr.ib(unit=u.Mpc/u.h)
    m = uttr.ib(unit=u.M_sun)

    _len = uttr.ib(init=False)

    def __attrs_post_init__(self):
        lengths = {
            len(e) for e in (self.x, self.y, self.z, self.vx, self.vy, self.vz)
        }

        if len(lengths) != 1:
            raise ValueError("Arrays should be of the same size")

        super().__setattr__("_len", lengths.pop())

    def __len__(self):
        return self._len

    def __repr__(self):
        cls_name = type(self).__name__
        length = len(self)
        return f"<{cls_name} size={length}>"


