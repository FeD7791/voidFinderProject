import os
import ctypes
import numpy as np
from .. import spherical_voids





def spherical_void_finder(box):
	
	#Import library
	path = os.getcwd()
	print(path)
	
	clibrary = ctypes.CDLL(os.path.join(path,'voidfindertk','spherical','vf_lib.so'),mode=ctypes.RTLD_GLOBAL) 
	

	#Create Pointer
	arr_pointer_x = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_y = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_z = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_vx = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_vy = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_vz = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_m = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])

	#Create stucts --> class
	class voids(ctypes.Structure):
	    _fields_ = [("n_voids",ctypes.c_int),
			("Rad", ctypes.c_float),
		        #("Rini", ctypes.c_float),
		        #("Ini", ctypes.c_float * 1),
		        ("Pos", ctypes.c_float * 3),
		        ("Vel", ctypes.c_float * 3),
		        ("Dtype", ctypes.c_float),
		        ("Delta", ctypes.c_float),
		        ("Poisson", ctypes.c_float),
		        #("Dist4", ctypes.c_float),
		        #("ToF", ctypes.c_bool),
		        ("Nran", ctypes.c_int)
	]

	class voidArray(ctypes.Structure):
	    _fields_ = [("voids", voids * 1)]


	#Declare Input Pointers
	clibrary.execute_void_finder.argtypes = [
	    arr_pointer_x,arr_pointer_y,arr_pointer_z,
	    arr_pointer_vx,arr_pointer_vy,arr_pointer_vz,arr_pointer_m,ctypes.c_int]

	#Declare OUTPUT args
	clibrary.execute_void_finder.restype = ctypes.POINTER(voidArray)

	


	#Return void array
	va = clibrary.execute_void_finder(box.x,box.y,box.z,box.vx,box.vy,box.vz,box.m,len(box))
	
	va_arr = np.ctypeslib.as_array(va, shape=(10,))
	
	n_voids = va_arr[0][0][0][0]-1
	
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
	for i in range(0,n_voids,1):
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
	
	sp_void = spherical_voids.SphericalVoids(
		rad = radius,
		x_void = x_coord,
		y_void = y_coord,
		z_void = z_coord,
		vel_x_void = vx_coord,
		vel_y_void = vy_coord,
		vel_z_void = vz_coord,
		delta = delta,
		dtype = dtype,
		poisson = poisson,
		nran = nran
	)
	
	
	return sp_void
	



