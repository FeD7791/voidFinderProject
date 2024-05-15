import math
from astropy import units as u
import numpy as np
import uttr
import os
import ctypes
from ..models import ModelABC
from ..analysis_tools import join_box_void
from ..tools import ctypes_input_params_builder, ctypes_output_builder, process_output_from_finder, preprocess_data_box

class SphericalVF(ModelABC):
    def __init__(self):
        pass

    def preprocess(self, databox, **kwargs):
        databox = preprocess_data_box(databox=databox, **kwargs)
        return databox

    def model_find(self, llbox,**kwargs):
        sp_void = spherical_void_finder(llbox.box,**kwargs)
        return {'voids':sp_void}
    
    def mk_vbox(self, voids,llbox):
        voids = voids['voids']
        # databox.box.__dict__.pop('_len')
        # voids.__dict__.pop('_void_len')
        #sparse_m = Classifier(**databox.box.__dict__, **voids.__dict__)._sparse_matrix(0.0)
        box_void_sparse = join_box_void(llbox.box, voids, tol=0.0)
        return box_void_sparse


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
    
    def __eq__(self, other):
        return all(
            [
                np.array_equal(self.rad, other.rad),
                np.array_equal(self.x_void, other.x_void),
                np.array_equal(self.y_void, other.y_void),
                np.array_equal(self.z_void, other.z_void),
                np.array_equal(self.vel_x_void, other.vel_x_void),
                np.array_equal(self.vel_y_void, other.vel_y_void),
                np.array_equal(self.vel_z_void, other.vel_z_void),
                np.array_equal(self.delta, other.delta),
                np.array_equal(self.dtype, other.dtype),
                np.array_equal(self.poisson, other.poisson),
                np.array_equal(self.nran, other.nran),
            ]
        )
    
    # def _slice(self,min:int,max:int):
    #     sp_voids = {
    #         'rad':self.rad[min:max],
    #         'x_void':self.x_void[min:max],
    #         'y_void':self.y_void[min:max],
    #         'z_void':self.z_void[min:max],
    #         'vel_x_void':self.vel_x_void[min:max],
    #         'vel_y_void':self.vel_y_void[min:max],
    #         'vel_z_void':self.vel_z_void[min:max],
    #         'delta':self.delta[min:max],
    #         'dtype':self.dtype[min:max],
    #         'poisson':self.poisson[min:max],
    #         # self.dist4,
    #         'nran':self.nran[min:max],

    #     }
    #     return SphericalVoids(**sp_voids)

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

def spherical_void_finder(box,**kwargs):
    """Uses an external C library (`vf_lib.so`) to identify spherical voids in
    a simulation box.

    Parameters
    ----------
    box : object
    An object containing the following attributes:
        * x (numpy.ndarray): x-coordinates of particles in the box (float64).
        * y (numpy.ndarray): y-coordinates of particles in the box (float64).
        * z (numpy.ndarray): z-coordinates of particles in the box (float64).
        * vx (numpy.ndarray): x-velocity of particles in the box (float64).
        * vy (numpy.ndarray): y-velocity of particles in the box (float64).
        * vz (numpy.ndarray): z-velocity of particles in the box (float64).
        * m (numpy.ndarray): mass of particles in the box (float64).

    Returns
    -------
    spherical_voids.SphericalVoids
        An instance of the SphericalVoids class containing the identified
        voids' properties.

    Notes
    -----
    This function relies on the external C library `vf_lib.so` located in
    the "voidfindertk/spherical" directory relative to the current working
    directory. The library is not loaded or included here.

    The function extracts void properties from the C library output and
    populates a SphericalVoids object. Refer to the SphericalVoids docstring
    for details on the returned object's attributes.
    """
    params = ctypes_input_params_builder(
        'spherical',
        BoxSize = math.ceil(np.max(box.x.value)),
        **kwargs
        ) #Build input params for input of finder
    
    # Import library
    path = os.path.dirname(os.path.realpath(__file__))
    
    clibrary = ctypes.CDLL(
        os.path.join(path,"spherical", "vf_lib.so"),
        mode=ctypes.RTLD_GLOBAL,
    ) # Path of the library

    arr_pointers = 7*[np.ctypeslib.ndpointer(
        dtype=np.float64, ndim=1, flags=["CONTIGUOUS"]
    )] #pointers for box input x,y,z,vx,vy,vz,m
    
    # Create stucts --> class for output
    output_class = ctypes_output_builder('spherical')
    
    # Declare Input Pointers
    clibrary.execute_void_finder.argtypes = arr_pointers + [ctypes.c_int, ctypes.POINTER(params['InputParams_class'])]
    
    # Declare OUTPUT args
    clibrary.execute_void_finder.restype = ctypes.POINTER(output_class['voidArray_class'])

    # Call and Return array of voids
    va = clibrary.execute_void_finder(
        box.x, box.y, box.z, 
        box.vx, box.vy, box.vz, box.m, 
        len(box), params['InputParams_class'](**params['params_dict'])
    )
    output = process_output_from_finder('spherical', array_of_voids=va)
    sp_void = SphericalVoids(**output)

    return sp_void


        

