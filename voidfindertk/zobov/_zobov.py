# modulos de python
import ctypes
import os
import pathlib
import shutil
import subprocess
import tempfile


from astropy import units as u

import h5py

import numpy as np

import pandas as pd

import re

import scipy

import uttr

from ..models import DataBox, ModelABC
from ..tools import join_box_void

class _Paths:

    CURRENT = pathlib.Path(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    )
    ZOBOV = CURRENT / "src"
    LOADER_SO = ZOBOV / "zobov_loader.so"


class ZobovVF(ModelABC):

    def __init__(self, path=None, workdir=None, dtype=np.float32):
        self._loader_so = pathlib.Path(
            _Paths.LOADER_SO if path is None else pathlib.Path(path)
        )

        self._workdir = pathlib.Path(
            tempfile.mkdtemp(prefix=f"vftk_{type(self).__name__}_")
            if workdir is None
            else workdir
        )

        self._dtype = dtype

    @property
    def path_zobov(self):
        return self._loader_so.parent

    def preprocess(self, databox):
        # Create input binary files for Zobov finder

        # Declare library path
        clibrary = ctypes.CDLL(str(self._loader_so), mode=ctypes.RTLD_GLOBAL)

        # Create Input Pointers for x,y,z,vx,vz,vy,m
        arr_pointer = 7 * [
            np.ctypeslib.ndpointer(
                dtype=np.float64, ndim=1, flags=["CONTIGUOUS"]
            )
        ]

        # Declare Input Pointers type
        clibrary.c_binary_writter.argtypes = arr_pointer + [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_char_p,
        ]
        # Fill Input
        clibrary.c_binary_writter(
            databox.x,
            databox.y,
            databox.z,
            databox.vx,
            databox.vy,
            databox.vz,
            databox.m,
            len(databox),
            os.path.join(self.path_zobov, "tracers_zobov.raw").encode("utf-8"),
            os.path.join(self.path_zobov, "tracers_zobov.txt").encode("utf-8"),
        )
        return DataBox(databox)

    def model_find(self, llbox):
        sp_void = zobov_void_finder(llbox.box)
        #sp_void._tracers_in_void = calculate_tracers_inside_void(
        #    llbox.box, sp_void
        #)
        #$return {"voids": sp_void}
        return 0
        
    def mk_vbox(self, voids, llbox):
        voids = voids["voids"]
        box_void_sparse = join_box_void(llbox.box, voids, tol=0.0)
        return box_void_sparse

    def get_void_mass(self, voids, llbox):
        mass_list = []
        for i in range(len(voids)):
            # voids._tracers_in_void[i] = array of indexes of tracers in void
            mass_list.append(sum(llbox.box.m.value[voids._tracers_in_void[i]]))
        return mass_list


def _float32_converter(v):
    return np.asarray(v, dtype=np.float32)


@uttr.s(repr=False)
class ZobovVoids:
    Void_number = uttr.ib(converter=np.array)
    File_void_number = uttr.ib(converter=np.array)
    CoreParticle = uttr.ib(converter=np.array)
    CoreDens = uttr.ib(converter=np.array)
    ZoneVol = uttr.ib(converter=np.array)
    Zone_number_part = uttr.ib(converter=np.array)
    Void_number_Zones = uttr.ib(converter=np.array)
    VoidVol = uttr.ib(converter=np.array, unit=u.Mpc**3)
    Void_number_Part = uttr.ib(converter=np.array)
    VoidDensContrast = uttr.ib(converter=np.array)
    VoidProb = uttr.ib(converter=np.array)
    _void_len = uttr.ib(init=False)
    _tracers_in_void = uttr.ib(
        init=False, default=None
    )  # Provide tracers in void

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        lengths = set()
        for e in (
            self.Void_number,  # Rank of the void, in decreasing order of VoidDensContrast.
            self.File_void_number,  # Number of the void, in previous files
            self.CoreParticle,  # Index --Related to Box-- of the core particle of the void
            self.CoreDens,  # Density of the core particle
            self.ZoneVol,  # Volume of the central zone of the void
            self.Zone_number_part,  # Number of particles in the central zone of the void
            self.Void_number_Zones,  # Number of zones in the void
            self.VoidVol,  # Volume of the void
            self.Void_number_Part,  # Number of particles in the void
            self.VoidDensContrast,  # Contrast density of the void
            self.VoidProb,  # The probability that that DensContrast would arise from Poisson noise
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
    output = pd.read_csv(filename, sep="\s+", skiprows=2)
    output.columns = [
        "Void_number",
        "File_void_number",
        "CoreParticle",
        "CoreDens",
        "ZoneVol",
        "Zone_number_part",
        "Void_number_Zones",
        "VoidVol",
        "Void_number_Part",
        "VoidDensContrast",
        "VoidProb",
    ]
    zobov = ZobovVoids(**output.to_dict(orient="list"))
    return zobov


def zobov_void_finder(box, **kwargs):
    #Params
    kwargs.setdefault('buffer_size', 0.08)
    kwargs.setdefault('box_size', box.size())
    kwargs.setdefault('number_of_divisions',2)
    kwargs.setdefault('delete_files', True)
    kwargs.setdefault('density_threshold', 0.2)

    #Paths
    path = _Paths.CURRENT / "src"/"src"
    


   
    #shutil.move(os.path.join(path_zobov,'tracers_zobov.raw'),path_src)
    #os.chdir(path)
    print("#######",Path.cwd())
    #Runing Zobov 
    subprocess.run(["./vozinit", "../tracers_zobov.raw", 
                    str(kwargs['buffer_size']), 
                    str(kwargs['box_size']),
                    str(kwargs['number_of_divisions']),
                    'output_vozinit'], cwd=str(path)) #Specify cwd as the working directory
    subprocess.run(["./scroutput_vozinit"], cwd=str(path))
    subprocess.run(["./jozov", 
                    "adjoutput_vozinit.dat", 
                    "voloutput_vozinit.dat", 
                    'out_particle_zone.dat',
                    'out_zones_in_void.dat',
                    'out_text_file.dat',
                    str(kwargs['density_threshold'])], cwd=str(path))

    #Delete unused files
    if kwargs['delete_files']: #provide delete_files as false to preserve files
        subprocess.Popen((
            'find', '.', 
            '-type', 'f', 
            '-name', 'part.output_vozinit*', 
            '-exec', 'rm','{}', ';')) #Delete part... files
        #Delete other binary unused files
        files_to_remove = [
            '../tracers_zobov.raw',
            "adjoutput_vozinit.dat",
            "voloutput_vozinit.dat",
            #'out_particle_zone.dat',
            #'out_zones_in_void.dat',
            'scroutput_vozinit']
        for f in files_to_remove:
            try:
                subprocess.Popen(["rm",f], cwd=str(path))
                #os.remove(f)
            except FileNotFoundError:
                print(f'File {f} not found')
    #Output Results
    zobov_voids = read_zobov_output(str(path/'out_text_file.dat'))   
    return zobov_voids 


def calculate_tracers_inside_void(box, voids, **kwargs):
    kwargs.setdefault("hdf5", True)
    xyz_tracers = np.array(
        [
            np.array([box.x.value[i], box.y.value[i], box.z.value[i]])
            for i in range(len(box))
        ]
    )
    # filter this array based on the ones that are centers according to voids
    void_centers = [
        xyz_tracers[i] for i in voids.CoreParticle
    ]  # Centers of the void

    # remove the void centers from xyz_tracers
    ##############################################################################
    # # Reshape filter_arr to match the dimensionality of a for comparison
    # filter_arr = void_centers.reshape(-1)  # Flatten filter_arr

    # # Use np.in1d to check for membership, ravel to flatten the result
    # not_in_filter = xyz_tracers.ravel()[~np.in1d(xyz_tracers.ravel(), filter_arr)]

    # # Reshape the result to maintain the original dimensionality
    # xyz_tracers = not_in_filter.reshape(-1, 3)
    ##############################################################################
    if kwargs["hdf5"]:
        print("ATTEMPTING TO CALCULATE TRACERS")
        path = os.path.dirname(os.path.realpath(__file__))
        file = h5py.File(
            os.path.join(path, "tracers_in_voids.h5"), "w"
        )  # save in zobovvf
        for i in range(len(void_centers)):
            # distance from center to each particle : row[i] = [dist(xyz_void(i),box_xyz(i))]
            d = scipy.spatial.distance.cdist([void_centers[i]], xyz_tracers)
            # asociate index to each particle distances[i] = (i ,dist(xyz_void(i),box_xyz(i)))
            distances = [list(enumerate(arr)) for arr in d]
            sorted_distances = [
                sorted(dist, key=lambda x: x[1]) for dist in distances
            ]  # sort the array ascending
            # get the number of particles in each void
            n_tracer_in_voids = voids.Void_number_Part
            # keep the lowest n_tracer_in_voids[i] from sorted_distances

            # sd = [sorted_distances[i][1: n_tracer_in_voids[i]+1] for i in range(len(sorted_distances))]
            sd = [
                sorted_distances[j][1 : n_tracer_in_voids[i] + 1]
                for j in range(len(sorted_distances))
            ]

            file.create_dataset(f"{i}", data=sd)

        # Get indexes, returns list of list of indexes
        file.close()
        file2 = h5py.File(os.path.join(path, "tracers_in_voids.h5"), "r")

        index = [
            list(list(zip(*value[:][0]))[0])
            for key, value in dict(file2).items()
        ]
        index = [
            list(np.array(arr, dtype=int)) for arr in index
        ]  # transform in array of integers

        file2.close()

    return index
    
def find_zobov_tracers():
    path = _Paths.CURRENT / "src"
    subprocess.run(['./lectura_zones'], cwd=str(path))
    subprocess.run(['./particle_zones'], cwd=str(path))
    with open('txt_out_particle_zone2.txt','r') as f:
        out = f.readlines()
    
    with open('txt_out_zones_in_void.txt','r') as f:
        out2 = f.readlines()
    e = [
    {
    'n_zone':np.int32(out[i].split(" ")[2].split())[0], 
    'particles':np.int32(out[i+2].split(" ")[:-1])} for i in range(len(out)) if re.match(r"^ zone", out[i])]
    
    f = [{'void':np.int32(out2[i].split(" ")[3:4])[0], 'zones':np.int32(out2[i].split(" ")[5:-1])} for i in range(3,len(out2))]
    
      
      
      
      
      
