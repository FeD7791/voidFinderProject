import os
import pathlib
import shutil
import tempfile

import attr
import numpy as np

import pandas as pd

from astropy import units as u

from ..utils import make_workdir

from . import _wrapper
from . import _postprocessing
from ..vfinder_abc import ModelABC


class _Paths:
    CURRENT = pathlib.Path(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    )
    # Path to the src folder of Zobov
    SVF = CURRENT / "popcorn_void_finder" / "Source"
    CONFFILE = SVF / "configuration"

class _FileNames:
    CONFIG = "vars.conf"
    TRSFILE = "trsfile.dat"
    SPHFILE = "sphfile.dat"
    POPFILE = "popfile.dat"
    RAWPOPFILE = "rawpopfile.dat"
    PAIRSFILE = "pairsfile.dat"

@attr.define
class PopCornVF(ModelABC):
    """
    PopCornVF class for void finding and analysis.

    Attributes
    ----------
    auxfiles : str
        Flag indicating the use of auxiliary files.
    boxsize : float
        Length of the box in the same units as the tracer input file.
    densth : float
        Density threshold for void identification.
    minradius : int
        Minimum radius of a sphere in input units.
    maxradius : int
        Maximum radius of a sphere in input units.
    massmin : float
        Minimum halo mass allowed (0 if not applicable).
    svf_path : pathlib.Path
        Path to the source directory of the SVF.
    mpi_flags : str or None
        MPI flags for parallel computation.
    workdir : pathlib.Path
        Path to the working directory.
    workdir_clean : bool
        Flag to clean the working directory on deletion.

    Methods
    -------
    __attrs_post_init__()
        Initializes paths for SVF and working directory.
    auxfiles
        Returns the flag for auxiliary files.
    boxsize
        Returns the length of the box.
    densth
        Returns the density threshold.
    minradius
        Returns the minimum radius for void detection.
    maxradius
        Returns the maximum radius for void detection.
    massmin
        Returns the minimum mass threshold.
    svf_path
        Returns the path to the SVF directory.
    mpi_flags
        Returns the MPI flags.
    workdir
        Returns the working directory path.
    workdir_clean
        Returns the flag for cleaning the working directory.
    _create_run_work_dir()
        Creates and returns a temporary working directory.
    __del__()
        Cleans up the temporary working directory if
        `workdir_clean` is True.
    preprocess(databox)
        Placeholder for preprocessing the DataBox object.
    model_find(databox)
        Finds voids using the provided DataBox object.
    build_voids(model_find_parameters)
        Builds voids and returns properties and tracers
        from the model find parameters.
    """
    _auxfiles = attr.field(default="true")# AUXILIARY FILES
    # INPUT PARAMETERS
    _boxsize = attr.field(default=1000.0)
    _densth = attr.field(default=-0.9)
    _minradius = attr.field(default=5)
    _maxradius = attr.field(default=100)
    _massmin = attr.field(default=0)
    # Path to Source folder
    _svf_path = attr.field(default=None)# Path to source directory of SVF
    # mpi flags:
    _mpi_flags = attr.field(default=None)
    # Path to working directory
    _workdir = attr.field(default=None)
    # Whether to clean or not the working directory
    _workdir_clean = attr.field(default=False)

    # Set path to extra files
    def __attrs_post_init__(self):

        # svf_path
        self._svf_path = pathlib.Path(
            _Paths.SVF if self._svf_path is None else self._svf_path
        )

        # Set workdir path
        self._workdir = pathlib.Path(
            tempfile.mkdtemp(prefix=f"svf_{type(self).__name__}_")
            if self._workdir is None
            else pathlib.Path(os.path.abspath(self._workdir))
        )


    # PROPERTIES ==============================================================




    @property
    def auxfiles(self):
        return self._auxfiles

    @property
    def boxsize(self):
        return self._boxsize

    @property
    def densth(self):
        return self._densth

    @property
    def minradius(self):
        return self._minradius

    @property
    def maxradius(self):
        return self._maxradius

    @property
    def massmin(self):
        return self._massmin

    @property
    def svf_path(self):
        return self._svf_path

    @property
    def mpi_flags(self):
        return self._mpi_flags

    @property
    def workdir(self):
        return self._workdir

    @property
    def workdir_clean(self):
        return self._workdir_clean

    # INTERNAL ================================================================
    # Directory Creator
    def _create_run_work_dir(self):
        run_workdir = make_workdir.create_run_work_dir(
            workdir_path=self._workdir
        )
        return run_workdir

    # Directory Cleaner
    def __del__(self):
        """
        Destructor that cleans up the temporary working directory
        if workdir_clean is True.
        """
        if self._workdir_clean:
            shutil.rmtree(self._workdir)

    def preprocess(self, databox):
        return databox

    def model_find(self, databox):
        """
        Runs the POPCORN void finder by creating the input file and the
        vars.conf parameter file and then runind the binary using command line
        instructions.

        Parameters
        ----------
            databox :  Object
            DataBox object that holds the Box object with tracers properties.

        Returns : Dict
        Dictionary with two parameters, 
            run_work_dir : path to the working directory.
            box : Box Object.

        Notes
        -----
        To run POPCORN void finder an input tracer and a configuration file
        are needed. The input tracer file is built using the Box Object.
        The configuration file is built using the parameters of the class.

        """
        # Retrieve box from DataBox object
        box = databox.box
        # create the sandbox
        run_work_dir = self._create_run_work_dir()
        # Create config file on Workdir
        _wrapper.config_file_maker(
            # Files
            trsfile=str(run_work_dir / _FileNames.TRSFILE),
            filefmt="ASCII",
            num_file=str(1),
            sphfile=str(run_work_dir / _FileNames.SPHFILE),
            popfile=str(run_work_dir / _FileNames.POPFILE),
            auxfiles=str(self._auxfiles),
            rawpopfile=str(run_work_dir / _FileNames.RAWPOPFILE),
            pairsfile=str(run_work_dir / _FileNames.PAIRSFILE),
            #Parameters
            boxsize=str(self._boxsize),
            densth=str(self._densth),
            minradius=str(self._minradius),
            maxradius=str(self._maxradius),
            massmin=str(self._massmin),
            eps=str(1e-5),
            path=str(run_work_dir / _FileNames.CONFIG),  # Workdir path
        )
        # Generate dataset file from box
        _wrapper.popcorn_svf_input_data_builder(
            box=box, file_path=str(run_work_dir / _FileNames.TRSFILE)
        )# Save File to workdir
        # Run Void Finder
        _wrapper.spherical_popcorn_void_finder(
            mpi_flags=self._mpi_flags,
            bin_path=_Paths.SVF,
            conf_file_path=run_work_dir,
            work_dir_path=run_work_dir,
        )
        return {"run_work_dir": run_work_dir, "box":box}

    def build_voids(self, model_find_parameters):
        """
        Postprocess the outputs of POPCORN (spherical) void finder to get the
        list of tracers inside each void (if any) and properties foud by this
        method.

        Parameters
        ----------
            model_find_parameters : Dict
            Parameters got from the model_find step.

        Returns
        -------
            tracers_in_voids : tuple of
            List of indexes of tracers (Regarding to the Box object index)
            inside each void.

            extra : Dict
            Dictionary with extra parameters, varying from properties to path
            directorys.
        """
        # Retrieve box from DataBox object
        box = model_find_parameters["box"]
        # Get current working directory
        run_work_dir = model_find_parameters["run_work_dir"]
        # Get void Properties
        properties = _postprocessing.get_void_properties(
            popcorn_output_file_path=str(run_work_dir / _FileNames.SPHFILE)
        )
        # Get tracers in void
        tracers_in_voids = _postprocessing.get_tracers_in_voids(
            box=box,
            popcorn_output_file_path=str(run_work_dir / _FileNames.SPHFILE)
        )
        # Build extra
        extra = {
            "svf_voids_properties":properties,
            "files_directory_path":run_work_dir
        }
        return tuple(tracers_in_voids), extra


# @uttr.s(repr=False)
# class PopCornSVF:
#     radius = uttr.ib(converter=np.array,unit=u.Mpc)
#     x_void = uttr.ib(converter=np.array,unit=u.Mpc)
#     y_void = uttr.ib(converter=np.array,unit=u.Mpc)
#     z_void = uttr.ib(converter=np.array,unit=u.Mpc)
#     delta = uttr.ib(converter=np.array)

# def popcorn_svf_input_data_builder(box, **kwargs):
#     kwargs.setdefault("sep",'\t')
#     kwargs.setdefault("index",False)
#     kwargs.setdefault("header",False)
#     path_file = os.path.dirname(os.path.realpath(__file__))
#     df = pd.DataFrame(box.__dict__)
#     df.drop(labels=['_len'], axis=1,inplace=True)
#     df.to_csv(
#         os.path.join(path_file,'popcorn','popcorn_void_finder','Source','box.dat'),
#         **kwargs) #Save the box in a file

# def config_file_maker(params, path):
#     '''
#     path: place where the input file is (always has vars.conf as name)
#     '''
#     params = dict([(key,str(value)) for key,value in params.items()]) #Transform non str items in strings
#     config_element = ConfigParser()
#     config_element.optionxform = str
#     config_element['INPUT_PARAMS'] = params
#     with open(os.path.join(path,'vars.conf'), 'w') as conf:
#         config_element.write(conf)

# def spherical_popcorn_void_finder():
#     #Paths
#     path = os.path.dirname(os.path.realpath(__file__))
#     path_source = os.path.join(path,'popcorn','popcorn_void_finder','Source')
#     os.chdir(path_source)
#     subprocess.run(["mpirun", "-np", "1","--bind-to","none", "./svf", "./configuration/vars.conf"])

# class SphericalPopCornVF(ModelABC):
#     def __init__(self):
#         pass

#     def preprocess(self, databox, **kwargs):
#         box = databox.box
#         popcorn_svf_input_data_builder(box, **kwargs)
#         return databox

#     def model_find(self, llbox, **kwargs):
#         #Paths
#         path = os.path.dirname(os.path.realpath(__file__))
#         path_config_file = os.path.join(
#             path,
#             'popcorn',
#             'popcorn_void_finder',
#             'Source',
#             'configuration')
#         #Default Params
#         kwargs.setdefault('TRSFILE',
#                           os.path.join(
#                               path,
#                               'popcorn',
#                               'popcorn_void_finder',
#                               'Source',
#                               'box.dat'))
#         kwargs.setdefault('FILEFMT','ASCII')
#         kwargs.setdefault('NUM_FILE', 1)
#         kwargs.setdefault('POPFILE', '../MWE/clean_popvds.dat')
#         kwargs.setdefault('SPHFILE', '../MWE/sphvds.dat')
#         kwargs.setdefault('AUXFILES' , 'true')
#         kwargs.setdefault('RAWPOPFILE','../MWE/raw_popvds.dat')
#         kwargs.setdefault('RAWSPHFILE', '../MWE/raw_sphv.dat')
#         kwargs.setdefault('PAIRSFILE', '../MWE/pairs_raw_popvds.dat')
#         kwargs.setdefault('BOXSIZE' , 1000.0)
#         kwargs.setdefault('DENSTH', -0.9)
#         kwargs.setdefault('MINRADIUS', 5)
#         kwargs.setdefault('MAXRADIUS', 100)
#         kwargs.setdefault('EPS', 0.00001)
#         kwargs.setdefault('MASSMIN', 50)

#         # Create config file
#         config_file_maker(params=kwargs, path=path_config_file)
#         # Use finder
#         spherical_popcorn_void_finder()
#         # Get output data
#         sphvds = pd.read_csv(kwargs['SPHFILE'], sep='\s+')
#         sphvds.columns = [
#             'index',
#             'radius',
#             'x_center',
#             'y_center',
#             'z_center',
#             'delta'
#             ] #preguntar a dante cual es la salida
#         # Tansform output to PopCornSVF class
#         sphvds.drop(labels=['index'], axis=1, inplace=True)
#         popcorn_svf = PopCornSVF(**sphvds.to_dict(orient='list'))
#         return popcorn_svf

#     def mk_vbox(self, databox,voids,sparse,llbox):
#         voids = voids['voids']
#         sparse_m = join_box_void(databox,voids,tol=0.0)
#         return sparse_m


# @uttr.s(repr=False)
# class PopCornVoids:
#     void_id = uttr.ib()
#     number_members = uttr.ib()
#     void_volume = uttr.ib(converter=np.array,unit=u.Mpc**3)
#     number_subhalos = uttr.ib()
#     internal_flag = uttr.ib()
#     spheres = uttr.ib()
#     ID_subhalo = uttr.ib()
#     _void_len = uttr.ib(init=False)

#     def __attrs_post_init__(self):
#         """Post init method.

#         Checks that the lenght of the inputs are the same
#         """
#         lengths = set(())
#         for e in (
#             self.void_id,
#             self.number_members,
#             self.void_volume ,
#             self.number_subhalos,
#             self.internal_flag,
#             self.spheres,
#             self.ID_subhalo
#         ):
#             lengths.add(len(e))

#         if len(lengths) != 1:
#             raise ValueError("Arrays should be of the same size")

#         super().__setattr__("_void_len", lengths.pop())

#     def __len__(self):
#         """Length method.

#         Returns
#         -------
#             int
#                 the number of elements in SphericalVoids
#         """
#         return self._void_len

#     def _slice(self,min,max):
#         popvoid = {
#             'void_id':self.void_id[min:max],
#             'number_members':self.number_members[min:max],
#             'void_volume':self.void_volume[min:max],
#             'number_subhalos':self.number_subhalos[min:max],
#             'internal_flag':self.internal_flag[min:max],
#             'spheres':self.spheres[min:max],
#             'ID_subhalo':self.ID_subhalo[min:max]
#         }
#         return PopCornVoids(**popvoid)

#     def __repr__(self):
#         """Representation method.

#         Returns
#         -------
#             str
#                 Name plus number of points in the box
#         """
#         cls_name = type(self).__name__
#         length = len(self)
#         return f"<{cls_name} size={length}>"
