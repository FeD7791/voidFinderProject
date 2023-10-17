import pandas as pd
import numpy as np
from tools import to_np_arr, slicer
from attrs import asdict, define, field, converters


@define
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

    x = field(converter=to_np_arr)
    y = field(converter=to_np_arr)
    z = field(converter=to_np_arr)
    vx = field(converter=to_np_arr)
    vy = field(converter=to_np_arr)
    vz = field(converter=to_np_arr)

    def build_df(self):
        data_df = np.vstack((self.x, self.y, self.z, self.vx, self.vy, self.vz))
        data_df = pd.DataFrame(data_df.transpose())
        data_df.columns = ["x", "y", "z", "vx", "vy", "vz"]
        return data_df

    def _check_length(self):
        if (
            len(self.x) == len(self.y)
            and len(self.x) == len(self.z)
            and len(self.vx) == len(self.x)
            and len(self.x) == len(self.vy)
            and len(self.vz) == len(self.x)
        ):
            pass
        else:
            print("Arrays should be of the same size")
            raise ValueError

    def slice(self, n, parameter):
        try:
            if parameter == "pos":
                columns = ["x", "z", "y", "distance"]
                slices = slicer(self.x, self.y, self.z, n, columns)
            elif parameter == "vel":
                columns = ["vx", "vz", "vy", "v modulus"]
                slices = slicer(self.vx, self.vy, self.vz, n, columns)
            return slices
        except ValueError:
            return "Invalid Parameter input. parameter = {pos,vel} "

    def __attrs_post_init__(self):
        self._check_length()
