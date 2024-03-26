from astropy import units as u
import numpy as np
import uttr
from abc import ABC
import os
import ctypes
import scipy
from attrs import define
from .box import Box
from .models import ModelABC


class SphericalVF(ModelABC):
    def __init__(self):
        pass

    def preprocess(self, databox):
        return databox

    def model_find(self, llbox):
        sp_void = spherical_void_finder(llbox.box)
        return {'voids':sp_void}
    
    def mk_vbox(self, databox,voids,llbox):
        voids = voids['voids']
        databox.box.__dict__.pop('_len')
        voids.__dict__.pop('_void_len')
        #sparse_m = Classifier(**databox.box.__dict__, **voids.__dict__)._sparse_matrix(0.0)
        sparse_m = Classifier(**databox.box.__dict__, **voids.__dict__)
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

def spherical_void_finder(box):
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
    # Import library
    path = os.path.dirname(os.path.realpath(__file__))
    

    clibrary = ctypes.CDLL(
        os.path.join(path,"spherical", "vf_lib.so"),
        mode=ctypes.RTLD_GLOBAL,
    )

    # Create Pointer
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

    # Create stucts --> class
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
    ]

    # Declare OUTPUT args
    clibrary.execute_void_finder.restype = ctypes.POINTER(voidArray)

    # Return void array
    va = clibrary.execute_void_finder(
        box.x, box.y, box.z, box.vx, box.vy, box.vz, box.m, len(box)
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

@define(repr=False)
class Classifier(Box, SphericalVoids):
    """Inherits from both Box and ShpericalVoids classes and provides a method
    to create a sparse matrix representing the relationship between spherical
    voids and tracer particles within a simulation box.

    Methods
    -------
    _sparse_matrix(self, tol=0.0)
        Generates a sparse CSR matrix indicating which voids are potentially
        influencing each tracer particle based on their relative positions.

    Parameters
    ----------
    self : Classifier object
        The instance of the class.
    tol : float, optional
        Tolerance level for considering a void to be influencing a tracer
        particle. Defaults to 0.0.

    Returns
    -------
    scipy.sparse.csr_matrix
        A sparse CSR matrix where each row represents a tracer particle and
        each column represents a void. The value at a specific row-column
        intersection is 1.0 if the corresponding void is potentially
        influencing the tracer particle (Tracer inside the void) based on
        the distance between their centers and the void's radius,
        otherwise 0.0.
    """

    def _sparse_matrix(self, tol=0.0):
        tolerance = tol * np.ones(self.rad.shape)
        rad = np.array(self.rad) + tolerance
        pos_void = np.array([self.x_void, self.y_void, self.z_void])
        pos_void = np.transpose(pos_void)
        arr = []
        for i in range(self._len):
            pos_box = np.array(
                [[self.x[i].value, self.y[i].value, self.z[i].value]]
            )
            d = scipy.spatial.distance.cdist(
                pos_box, pos_void, metric="euclidean"
            )  # Calculates the distance from a tracer to each void center
            rad_center_distance_comparing = np.greater(rad, d[0]).astype(
                float
            )  # Compares radius of a void and distance from the void center to the particle 1 if rad>d # noqa: E501
            arr.append(rad_center_distance_comparing)
        # Transform to sparse matrix
        output_sparse_matrix = scipy.sparse.csr_matrix(arr)
        return output_sparse_matrix
    
    def _get_voids(self):
            
        
        return {
            'rad':self.rad,
            'x_void':self.x_void,
            'y_void':self.y_void,
            'z_void':self.z_void,
            'vel_x_void':self.vel_x_void,
            'vel_y_void':self.vel_y_void,
            'vel_z_void':self.vel_z_void,
            'delta':self.delta,
            'dtype':self.dtype,
            'poisson':self.poisson,
            'nran':self.nran
            
            }
            
    def _get_box(self):
        return{
        'x':self.x,
        'y':self.y,
        'z':self.z,
        'vx':self.vx,
        'vy':self.vy,
        'vz':self.vz,
        'm':self.m
        }

    def __repr__(self):
        """Representation method.

        Returns
        -------
            str
                Size of Box and calculated Voids
        """
        cls_name = type(self).__name__
        # length_box = self._len
        # length_voids = self._void_len
        return f"<{cls_name} >"
        

