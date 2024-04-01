import uttr
import numpy as np
from astropy import units as u
import os
import ctypes
import subprocess
import os
import shutil



@uttr.s(repr=False)
class ZobovVoids:
    Void_number = uttr.ib(converter=np.array)
    File_void_number = uttr.ib(converter=np.array)
    CoreParticle = uttr.ib(converter=np.array)
    CoreDens = uttr.ib(converter=np.array)
    ZoneVol = uttr.ib(converter=np.array)
    Zone_number_part = uttr.ib(converter=np.array)
    Void_number_Zones = uttr.ib(converter=np.array)
    VoidVol = uttr.ib(converter=np.array,unit=u.Mpc**3)
    Void_number_Part = uttr.ib(converter=np.array)
    VoidDensContrast = uttr.ib(converter=np.array)
    VoidProb = uttr.ib(converter=np.array)
    _void_len = uttr.ib(init=False)

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        lengths = set(())
        for e in (
            self.Void_number,
            self.File_void_number,
            self.CoreParticle,
            self.CoreDens,
            self.ZoneVol,
            self.Zone_number_part,
            self.Void_number_Zones,
            self.VoidVol,
            self.Void_number_Part,
            self.VoidDensContrast,
            self.VoidProb
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
    
    def _slice(self,min,max):
        zobovvoid = {
            'Void_number':self.Void_number[min:max],
            'File_void_number':self.File_void_number[min:max],
            'CoreParticle':self.CoreParticle[min:max],
            'CoreDens':self.CoreDens[min:max],
            'ZoneVol':self.ZoneVol[min:max],
            'Zone_number_part':self.Zone_number_part[min:max],
            'Void_number_Zones':self.Void_number_Zones[min:max],
            'VoidVol':self.VoidVol[min:max],
            'Void_number_Part':self.Void_number_Part[min:max],
            'VoidDensContrast':self.VoidDensContrast[min:max],
            'VoidProb':self.VoidProb[min:max]
        }
        return ZobovVoids(**zobovvoid)
    
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
    
def zobov_void_finder(box):
	path = os.path.dirname(os.path.realpath(__file__))
	path_zobov = os.path.join(path,'zobov')
	
	#Declaramos la ubicacion de la libreria
	clibrary = ctypes.CDLL(os.path.join(path,'zobov','zobov_loader.so') ,mode=ctypes.RTLD_GLOBAL) 

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
	    arr_pointer_vx,arr_pointer_vy,arr_pointer_vz,arr_pointer_m,ctypes.c_int,
        ctypes.c_char_p, ctypes.c_char_p]
	#Use function c_binary_writter
	clibrary.c_binary_writter(box.x,box.y,box.z,box.vx,box.vy,box.vz,box.m,len(box),
                os.path.join(path_zobov,'tracers_zobov.raw').encode('utf-8'),
                os.path.join(path_zobov,'tracers_zobov.txt').encode('utf-8')
                           )




	
	
	path_src = os.path.join(path_zobov,'src')
	shutil.move(os.path.join(path_zobov,'tracers_zobov.raw'),path_src)
	os.chdir(path_src)
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
               

