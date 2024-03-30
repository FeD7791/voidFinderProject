import math
from astropy import units as u
import numpy as np
import uttr
import os
import ctypes
from .models import ModelABC
from .analysis_tools import join_box_void

class SphericalVF(ModelABC):
    def __init__(self):
        pass

    def preprocess(self, databox):
        return databox

    def model_find(self, llbox):
        sp_void = spherical_void_finder(llbox.box)
        return {'voids':sp_void}
    
    def mk_vbox(self, databox,voids,sparse,llbox):
        voids = voids['voids']
        # databox.box.__dict__.pop('_len')
        # voids.__dict__.pop('_void_len')
        #sparse_m = Classifier(**databox.box.__dict__, **voids.__dict__)._sparse_matrix(0.0)
        sparse_m = join_box_void(databox,voids,tol=0.0)
        return sparse_m


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
    
    def _slice(self,min:int,max:int):
        sp_voids = {
            'rad':self.rad[min:max],
            'x_void':self.x_void[min:max],
            'y_void':self.y_void[min:max],
            'z_void':self.z_void[min:max],
            'vel_x_void':self.vel_x_void[min:max],
            'vel_y_void':self.vel_y_void[min:max],
            'vel_z_void':self.vel_z_void[min:max],
            'delta':self.delta[min:max],
            'dtype':self.dtype[min:max],
            'poisson':self.poisson[min:max],
            # self.dist4,
            'nran':self.nran[min:max],

        }
        return SphericalVoids(**sp_voids)

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

    # Input params
    params = {
    #Double Values
    "RadIncrement" : 0. ,
    "BoxSize": math.ceil(np.max(box.x.value)) ,
    "MaxRadiusSearch": 40.0,
    "ProxyGridSize": 5.0,
    "FracRadius" : 0.5 ,
    "DeltaThreshold": -0.9,  
    "DeltaSeed": -0.7,
    "OverlapTol":0,
    "Redshift" : 0.99 ,
    "OmegaMatter" : 0.25,
    "OmegaLambda" : 0.75,
    "Hubble" : 0.73,
    "FidOmegaMatter ": 0.2  ,
    "FidOmegaLambda" : 0.8 ,
    "FidHubble" : 0.7, 
    "MinProfileDist ": 0.5,
    "MaxProfileDist" : 3.0,
    "ScalePos":1,
    "ScaleVel" : 1,
    "InnerShellVel": 0.8,
    "OuterShellVel": 1.2,
    #Int Values
    "FormatTracers":0,
    "NumFiles":32,
    "NumRanWalk" : 75 , 
    "OMPcores":8,
    "RSDist" : 0,
    "GDist ": 0,
    "WriteProfiles" : 0 ,
    "NumProfileBins" : 100,   
         
    }
    for key,value in kwargs.items():
        params[key] = value

    # Import library
    path = os.path.dirname(os.path.realpath(__file__))
    

    clibrary = ctypes.CDLL(
        os.path.join(path,"spherical", "vf_lib.so"),
        mode=ctypes.RTLD_GLOBAL,
    )

    # Load Input Params
    params1 = dict(list(params.items())[0:21]) # Params that are double type
    params2 = dict(list(params.items())[21:29]) # Params that are int type

    # Build ctypes structure for params
    p1 = [(key,ctypes.c_double) for key,value in params1.items()] 
    p2 = [(key,ctypes.c_int) for key,value in params2.items()]

    class InputParams(ctypes.Structure):
        _fields_ = p1 + p2

    # Create Pointer for input box elements
    arr_pointer_x = np.ctypeslib.ndpointer(
        dtype=np.float64, ndim=1, flags=["CONTIGUOUS"]
    )
    arr_pointer_y = np.ctypeslib.ndpointer(
        dtype=np.float64, ndim=1, flags=["CONTIGUOUS"]
    )
    arr_pointer_z = np.ctypeslib.ndpointer(
        dtype=np.float64, ndim=1, flags=["CONTIGUOUS"]
    )
    arr_pointer_vx = np.ctypeslib.ndpointer(
        dtype=np.float64, ndim=1, flags=["CONTIGUOUS"]
    )
    arr_pointer_vy = np.ctypeslib.ndpointer(
        dtype=np.float64, ndim=1, flags=["CONTIGUOUS"]
    )
    arr_pointer_vz = np.ctypeslib.ndpointer(
        dtype=np.float64, ndim=1, flags=["CONTIGUOUS"]
    )
    arr_pointer_m = np.ctypeslib.ndpointer(
        dtype=np.float64, ndim=1, flags=["CONTIGUOUS"]
    )

    # Create stucts --> class for output
    class voids(ctypes.Structure):
        _fields_ = [
            ("n_voids", ctypes.c_int),
            ("Rad", ctypes.c_float),
            # ("Rini", ctypes.c_float),
            # ("Ini", ctypes.c_float * 1),
            ("Pos", ctypes.c_float * 3),
            ("Vel", ctypes.c_float * 3),
            ("Dtype", ctypes.c_float),
            ("Delta", ctypes.c_float),
            ("Poisson", ctypes.c_float),
            # ("Dist4", ctypes.c_float),
            # ("ToF", ctypes.c_bool),
            ("Nran", ctypes.c_int),
        ]

    class voidArray(ctypes.Structure):
        _fields_ = [("voids", voids * 1)]

    # Declare Input Pointers
    clibrary.execute_void_finder.argtypes = [
        arr_pointer_x,
        arr_pointer_y,
        arr_pointer_z,
        arr_pointer_vx,
        arr_pointer_vy,
        arr_pointer_vz,
        arr_pointer_m,
        ctypes.c_int,
        ctypes.POINTER(InputParams) #Input params
        
    ]

    # Declare OUTPUT args
    clibrary.execute_void_finder.restype = ctypes.POINTER(voidArray)

    # Call and Return array of voids
    va = clibrary.execute_void_finder(
        box.x, box.y, box.z, 
        box.vx, box.vy, box.vz, box.m, 
        len(box), InputParams(**params)
    )

    va_arr = np.ctypeslib.as_array(va, shape=(10,))

    n_voids = va_arr[0][0][0][0] - 1

    full_arr = np.ctypeslib.as_array(va, shape=(n_voids,))

    radius = []
    x_coord = []
    y_coord = []
    z_coord = []
    vx_coord = []
    vy_coord = []
    vz_coord = []
    dtype = []
    delta = []
    poisson = []
    nran = []
    for i in range(0, n_voids, 1):
        radius.append(full_arr[i][0][0][1])
        x_coord.append(full_arr[i][0][0][2][0])
        y_coord.append(full_arr[i][0][0][2][1])
        z_coord.append(full_arr[i][0][0][2][2])
        vx_coord.append(full_arr[i][0][0][3][0])
        vy_coord.append(full_arr[i][0][0][3][1])
        vz_coord.append(full_arr[i][0][0][3][2])
        dtype.append(full_arr[i][0][0][4])
        delta.append(full_arr[i][0][0][5])
        poisson.append(full_arr[i][0][0][6])
        nran.append(full_arr[i][0][0][7])

    sp_void = SphericalVoids(
        rad=radius,
        x_void=x_coord,
        y_void=y_coord,
        z_void=z_coord,
        vel_x_void=vx_coord,
        vel_y_void=vy_coord,
        vel_z_void=vz_coord,
        delta=delta,
        dtype=dtype,
        poisson=poisson,
        nran=nran,
    )

    return sp_void


        

