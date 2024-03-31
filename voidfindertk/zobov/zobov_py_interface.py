import os
import ctypes
import numpy as np
import subprocess
import os
import shutil
from tools import io

def zobov_void_finder(box):
	path = os.getcwd()
	
	
	#Declaramos la ubicacion de la libreria
	clibrary = ctypes.CDLL(os.path.join(path,'zobov','zobov_loader'),mode=ctypes.RTLD_GLOBAL) 

	#Create Pointer
	arr_pointer_x = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_y = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_z = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_vx = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_vy = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_vz = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])
	arr_pointer_m = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['CONTIGUOUS'])

	#Declare Input Pointers
	clibrary.c_binary_writter.argtypes = [
	    arr_pointer_x,arr_pointer_y,arr_pointer_z,
	    arr_pointer_vx,arr_pointer_vy,arr_pointer_vz,arr_pointer_m,ctypes.c_int]
	#Use function c_binary_writter
	clibrary.c_binary_writter(box.x,box.y,box.z,box.vx,box.vy,box.vz,box.m,len(box))




	
	
	path2 = os.path.join(path,'zobov','src')
	shutil.move(os.path.join(path,'tracers_zobov.raw'),path2)
	os.chdir(path2)
	#subprocess.Popen(["./vozinit", "tracers_zobov.raw", '0.08', '500','2','output_vozinit'], stdin=subprocess.PIPE)
	subprocess.run(["./vozinit", "tracers_zobov.raw", '0.08', '500','2','output_vozinit'])
	subprocess.run(["./scroutput_vozinit"])
	subprocess.run(["./jozov", "adjoutput_vozinit.dat", "voloutput_vozinit.dat", 'out_particle_zone.dat','out_zones_in_void.dat','out_text_file.dat','0.2'])



	#Delete unused files
	subprocess.Popen(('find', '.', '-type', 'f', '-name', 'part.output_vozinit*', '-exec', 'rm','{}', ';')) #Delete part... files
	#Delete other binary unused files

	files_to_remove = ['tracers_zobov.raw',"adjoutput_vozinit.dat","voloutput_vozinit.dat",'out_particle_zone.dat','out_zones_in_void.dat','scroutput_vozinit']
	for f in files_to_remove:
		try:
			os.remove(f)
		except FileNotFoundError:
			print(f'File {f} not found')
	#Output Results
	
