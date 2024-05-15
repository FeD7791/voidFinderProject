import uttr
import numpy as np
from astropy import units as u
import pandas as pd
import os
import ctypes
import subprocess
import shutil
from ..models import ModelABC
from ..analysis_tools import join_box_void

class ZobovVF(ModelABC):
    def __init__(self):
        pass

    def preprocess(self, databox):
        '''
            This part will create a couple of files: 
            tracers_zobov.raw (binary file) and 
            tracers_zobov.dat a csv file
            tracers_zobov.raw is needed to feed zobov's
            finder input and make it work. To create the files
            a c script is used , it just takes the box reads it then
            saves it in the binary file. This method will return databox
            (which contains the box) as a dummy return (cause it will not be
            used for anything)
        '''
        #Create input binary files for Zobov finder
        path = os.path.dirname(os.path.realpath(__file__))
        path_zobov = os.path.join(path,'zobov')

        #Declare library path
        clibrary = ctypes.CDLL(
            os.path.join(path,'zobov','zobov_loader.so') ,
            mode=ctypes.RTLD_GLOBAL
            ) 

        #Create Input Pointers for x,y,z,vx,vz,vy,m
        arr_pointer = 7*[
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, 
                                   flags=['CONTIGUOUS'])]

        #Declare Input Pointers type
        clibrary.c_binary_writter.argtypes = arr_pointer + [
            ctypes.c_int,ctypes.c_char_p, ctypes.c_char_p]
        #Fill Input
        clibrary.c_binary_writter(
            databox.box.x,
            databox.box.y,
            databox.box.z,
            databox.box.vx,
            databox.box.vy,
            databox.box.vz,
            databox.box.m,
            len(databox.box),
            os.path.join(path_zobov,'tracers_zobov.raw').encode('utf-8'),
            os.path.join(path_zobov,'tracers_zobov.txt').encode('utf-8')
                            )
        return databox

    def model_find(self, llbox):
        sp_void = zobov_void_finder(llbox.box)
        return {'voids':sp_void}
    
    def mk_vbox(self, voids,llbox):
        voids = voids['voids']
        # databox.box.__dict__.pop('_len')
        # voids.__dict__.pop('_void_len')
        #sparse_m = Classifier(**databox.box.__dict__, **voids.__dict__)._sparse_matrix(0.0)
        box_void_sparse = join_box_void(llbox.box, voids, tol=0.0)
        return box_void_sparse

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
            self.Void_number, #Rank of the void, in decreasing order of VoidDensContrast.
            self.File_void_number, #Number of the void, in previous files
            self.CoreParticle, # Index --Related to Box-- of the core particle of the void 
            self.CoreDens, # Density of the core particle
            self.ZoneVol, # Volume of the central zone of the void
            self.Zone_number_part, # Number of particles in the central zone of the void
            self.Void_number_Zones, # Number of zones in the void
            self.VoidVol, # Volume of the void
            self.Void_number_Part, # Number of particles in the void
            self.VoidDensContrast, # Contrast density of the void
            self.VoidProb # The probability that that DensContrast would arise from Poisson noise
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
    
    # def _slice(self,min,max):
    #     zobovvoid = {
    #         'Void_number':self.Void_number[min:max],
    #         'File_void_number':self.File_void_number[min:max],
    #         'CoreParticle':self.CoreParticle[min:max],
    #         'CoreDens':self.CoreDens[min:max],
    #         'ZoneVol':self.ZoneVol[min:max],
    #         'Zone_number_part':self.Zone_number_part[min:max],
    #         'Void_number_Zones':self.Void_number_Zones[min:max],
    #         'VoidVol':self.VoidVol[min:max],
    #         'Void_number_Part':self.Void_number_Part[min:max],
    #         'VoidDensContrast':self.VoidDensContrast[min:max],
    #         'VoidProb':self.VoidProb[min:max]
    #     }
    #     return ZobovVoids(**zobovvoid)
    
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
    
def read_zobov_output(filename):
    output = pd.read_csv(filename, sep='\s+', skiprows=2)
    output.columns = [
        'Void_number',
        'File_void_number',
        'CoreParticle',
        'CoreDens',
        'ZoneVol',
        'Zone_number_part',
        'Void_number_Zones',
        'VoidVol',
        'Void_number_Part',
        'VoidDensContrast',
        'VoidProb'
    ]
    zobov = ZobovVoids(**output.to_dict(orient='list'))
    return zobov
    


def zobov_void_finder(box,**kwargs):
    #Params
    kwargs.setdefault('buffer_size', 0.08)
    kwargs.setdefault('box_size', box.size())
    kwargs.setdefault('number_of_divisions',2)
    kwargs.setdefault('delete_files', True)
    kwargs.setdefault('density_threshold', 0.2)

    #Paths
    path = os.path.dirname(os.path.realpath(__file__))
    path_zobov = os.path.join(path,'zobov')
    


    path_src = os.path.join(path_zobov,'src')
    shutil.move(os.path.join(path_zobov,'tracers_zobov.raw'),path_src)
    os.chdir(path_src)

    #Runing Zobov 
    subprocess.run(["./vozinit", "tracers_zobov.raw", 
                    str(kwargs['buffer_size']), 
                    str(kwargs['box_size']),
                    str(kwargs['number_of_divisions']),
                    'output_vozinit'])
    subprocess.run(["./scroutput_vozinit"])
    subprocess.run(["./jozov", 
                    "adjoutput_vozinit.dat", 
                    "voloutput_vozinit.dat", 
                    'out_particle_zone.dat',
                    'out_zones_in_void.dat',
                    'out_text_file.dat',
                    str(kwargs['density_threshold'])])

    #Delete unused files
    if kwargs['delete_files']: #provide delete_files as false to preserve files
        subprocess.Popen((
            'find', '.', 
            '-type', 'f', 
            '-name', 'part.output_vozinit*', 
            '-exec', 'rm','{}', ';')) #Delete part... files
        #Delete other binary unused files
        files_to_remove = [
            'tracers_zobov.raw',
            "adjoutput_vozinit.dat",
            "voloutput_vozinit.dat",
            'out_particle_zone.dat',
            'out_zones_in_void.dat',
            'scroutput_vozinit']
        for f in files_to_remove:
            try:
                os.remove(f)
            except FileNotFoundError:
                print(f'File {f} not found')
    #Output Results
    zobov_voids = read_zobov_output(os.path.join(path_src,'out_text_file.dat'))   
    return zobov_voids     

def calculate_tracers_inside_void(box,voids):