"""
Module that holds functions and methods that are used to run the ZOBOV python
wrapper methods in a coherent step by step.

    ZobovVF is the main Model object that represent a particular run of ZOBOV
    void finder. It mainly consists of two steps:

    a. preprocess step used to preprocess any input data needed for the
    algorythm. For this method none pre process is needed
    b. model_find: Consist of all steps needed to run zobov and the python
    wrappers used to represent them. For reference of the steps see ZobovVF
    Class in this module.
"""

import datetime as dt
import os
import pathlib
import shutil
import tempfile

import numpy as np

from . import _postprocessing
from . import _wrapper as _wrap
from ..vfinder_abc import ModelABC


class _Names:
    OUTPUT_VOZINIT = "output_vozinit"
    OUTPUT_JOZOV_VOIDS = "output_txt"
    PARTICLES_IN_ZONES = "part_vs_zone"
    ZONES_IN_VOID = "zones_vs_voids"


class _Files:
    """
    Name of the Box parsed to raw files. This files always are parsed and
    saved in files with the same names (see write_input module in _wrapper)
    """

    TRACERS_RAW = "tracers_zobov.raw"
    TRACERS_TXT = "tracers_zobov.txt"
    PARTICLES_VS_ZONES_RAW = f"{_Names.PARTICLES_IN_ZONES}.dat"
    PARTICLES_VS_ZONES_ASCII = f"{_Names.PARTICLES_IN_ZONES}_ascii.txt"
    OUTPUT_JOZOV_VOIDS_DAT = f"{_Names.OUTPUT_JOZOV_VOIDS}.dat"
    ZONES_VS_VOID_RAW = f"{_Names.ZONES_IN_VOID}.dat"
    ZONES_VS_VOID_ASCII = f"{_Names.ZONES_IN_VOID}_ascii.txt"


class _ExecutableNames:
    ZOBOV_LOADER_BIN = "zobov_loader.so"
    TRACERS_IN_ZONES_BIN = "tracers_in_zones.so"
    ZONES_IN_VOIDS_BIN = "zones_in_void.so"


class _Paths:
    """
    Class that holds paths of reference to the current file and
    ZOBOV's src directory
    """

    CURRENT = pathlib.Path(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    )
    ZOBOV = CURRENT / "src"  # Path to the src folder of Zobov


class ZobovVF(ModelABC):
    """
    ZobovVF class for running ZOBOV Void Finder.

    This class provides methods to preprocess data and execute the ZOBOV
    Void Finder algorithm on a given data box.

    Parameters
    ----------
    buffer_size : float, optional
        Buffer size for ZOBOV (default is 0.08). The buffer size sets the size
        in units such that the box size of the data cube is 1, of the buffer
        around each sub-box when calculating the Voronoi diagram.
    box_size : int, optional
        Range of positions of particles in each dimension (default is 500).
    number_of_divisions : int, optional
        Number of divisions in each dimension of the box (default is 2).
    density_threshold : int, optional
        The density threshold is an optional parameter, which can limit the
        growth of voids into high-density regions. (default is 0).
    zobov_path : str or None, optional
        Path to ZOBOV executable (default is None, uses internal path).
    workdir : str or None, optional
        Temporary working directory path (default is None, creates a new temp
        directory).
    workdir_clean : bool, optional
        Whether to clean up the working directory on deletion (default is
        False).
    dtype : numpy.dtype, optional
        Data type used for computations (default is np.float32).

    Attributes
    ----------
        _buffer_size : float
            Input parameter of vozinit in ZOBOV void finder:
            The buffer size sets the size, in units such that the box size of
            the data cube is 1, of the buffer around each sub-box when
            calculating the Voronoi diagram.
        _box_size : float
            Input parameter of vozinit in ZOBOV void finder:
            The range of positions of particles in each dimension
        _number_of_divisions : int
            Input parameter of vozinit in ZOBOV void finder:
            (default 2) -- the no. of partitions in each dimension; must be at
            least 2 (giving 8 sub-boxes)
        _density_threshold : float
            Input parameter of vozinit in ZOBOV void finder:
        _zobov_path : Pathlib.path
        _workdir : Pathlib.path
        _workdir_clean : bool
        _dtype : numpy.dtype


    Methods
    -------
    preprocess(databox)
        Placeholder method for data preprocessing.
    model_find(databox)
        Executes the ZOBOV Void Finder algorithm on the provided DataBox
        object.
        This step follows these steps:
            1. Build the input data from the input box. This step will parse
            the box data to a raw file that the next step needs
            2. Run ZOBOV's vozinit executable using the input params and the
            tracers input file build in the last step. As the process ends
            an script file will be created (see run_vozinit).
            3. In this step the mentioned script will be run. It will result
            in the output of volume and adjacency files (see run_voztie).
            4. This step will run ZOBOV's jozov executable (see run_jozov)
            This step will result in the output of three files:
                - part_vs_zone.dat : Raw File containing the particles inside
                zones (see run_jozov)
                - zones_vsvoids.dat : Raw File containing the zones inside
                voids (see run_jozov)
                - output_txt.dat : Ascii File containing the voids properties
                (see run_jozov)
            5. This step will create the object Voids that contains the voids
            foud by the method and their properties.
    Notes
    -----
    The ZOBOV Void Finder is executed in several steps including VOZINIT,
    VOZSTEP, and JOZOV.
    """

    def __init__(
        self,
        *,
        buffer_size=0.08,
        box_size=500,
        number_of_divisions=2,
        density_threshold=0,
        zobov_path=None,
        workdir=None,
        workdir_clean=False,
        dtype=np.float32,
    ):

        self._buffer_size = buffer_size
        self._box_size = box_size
        self._number_of_divisions = number_of_divisions
        self._density_threshold = density_threshold

        self._zobov_path = pathlib.Path(
            _Paths.ZOBOV if zobov_path is None else zobov_path
        )
        # Create a workdir path to run ZOBOV
        self._workdir = pathlib.Path(
            tempfile.mkdtemp(prefix=f"vftk_{type(self).__name__}_")
            if workdir is None
            else pathlib.Path(os.path.abspath(workdir))
        )
        self._workdir_clean = bool(workdir_clean)

        self._dtype = dtype

    # PROPERTIES ==============================================================
    @property
    def buffer_size(self):
        return self._buffer_size

    @property
    def box_size(self):
        return self._box_size

    @property
    def number_of_divisions(self):
        return self._number_of_divisions

    @property
    def ensity_threshold(self):
        return self._density_threshold

    @property
    def zobov_path(self):
        return self._zobov_path

    @property
    def workdir(self):
        return self._workdir

    @property
    def workdir_clean(self):
        return self._workdir_clean

    @property
    def dtype(self):
        return self._dtype

    # INTERNAL ================================================================

    def _create_run_work_dir(self):
        """
        This method will create a temporal directory inside the working
        directory of the ZobovVF class workdir.

        Returns
        -------
            run_work_dir: pathlib.Path
                path of the work directoty
        """
        timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
        run_work_dir = pathlib.Path(
            tempfile.mkdtemp(suffix=timestamp, dir=self.workdir)
        )
        return run_work_dir

    def __del__(self):
        """
        Destructor that cleans up the temporary working directory
        if workdir_clean is True.
        """
        if self._workdir_clean:
            shutil.rmtree(self._workdir)

    def preprocess(self, databox):
        """
        Placeholder method for data preprocessing.

        Parameters
        ----------
        databox : object
            DataBox object containing data to be preprocessed.

        Returns
        -------
        object
            Preprocessed data.
        """
        return databox

    def model_find(self, databox):
        """
        Execute the ZOBOV Void Finder algorithm on the provided DataBox
        object.

        Parameters
        ----------
        databox : object
            DataBox object containing the data box to be analyzed.
        """
        # Retrieve box from DataBox object
        box = databox.box

        # create the sandbox
        run_work_dir = self._create_run_work_dir()

        # the tracers files
        tracers_raw_file_path = run_work_dir / _Files.TRACERS_RAW
        tracers_txt_file_path = run_work_dir / _Files.TRACERS_TXT

        # write the box in the files
        _wrap.write_input(
            box=box,
            path_executable=self._zobov_path
            / _ExecutableNames.ZOBOV_LOADER_BIN,
            raw_file_path=tracers_raw_file_path,
            txt_file_path=tracers_txt_file_path,
        )

        # VOZINIT =============================================================

        _wrap.run_vozinit(
            vozinit_dir_path=self._zobov_path / "src",
            input_file_path=tracers_raw_file_path,
            buffer_size=self.buffer_size,
            box_size=self.box_size,
            number_of_divisions=self.number_of_divisions,
            executable_name=_Names.OUTPUT_VOZINIT,
            work_dir_path=run_work_dir,
        )

        # VOZSTEP =============================================================
        # This step is mandatory if VOZINIT was run before

        _wrap.run_voz_step(
            preprocess_dir_path=run_work_dir,
            executable_name=_Names.OUTPUT_VOZINIT,
            work_dir_path=run_work_dir,
            voz_executables_path=_Paths.ZOBOV
            / "src",  # this is the path where voz1b1 and voztie exe are
        )

        # JOZOV ===============================================================
        _wrap.run_jozov(
            jozov_dir_path=_Paths.ZOBOV / "src",
            executable_name=_Names.OUTPUT_VOZINIT,
            output_name_particles_in_zones=_Names.PARTICLES_IN_ZONES,
            output_name_zones_in_void=_Names.ZONES_IN_VOID,
            output_name_text_file=_Names.OUTPUT_JOZOV_VOIDS,
            density_threshold=0,
            work_dir_path=run_work_dir,
        )
        return {"run_work_dir": run_work_dir}

    def build_voids(self, model_find_parameters):
        """
        This methods is used to build the final object Voids (see Voids class
        in this module). Each step will specify a mandatory attribute or
        method of Voids class.

        Parameters
        ----------
            model_find_parameters: Dictionary
                The dictionary holds some relevant properties from the
                model_find method (See model_find method in whithin this class
                ). Those properties are needed to run this module.
                These properties are:
                - run_work_dir: Directory path where the current run is
                performed.
                - databox : Object class DataBox (See DataBox in box module)
                with the tracers information.
        Returns
        -------
            voids : Voids class
                Object Voids that contains the voids found in this run
                alongside with their properties.
        """
        # Get current working directory
        run_work_dir = model_find_parameters["run_work_dir"]

        # Process 1:
        # a) Parse tracers in zones raw file in the work directory
        _postprocessing.parse_tracers_in_zones_output(
            executable_path=_Paths.ZOBOV
            / _ExecutableNames.TRACERS_IN_ZONES_BIN,
            input_file_path=run_work_dir
            / _Files.PARTICLES_VS_ZONES_RAW,
            output_file_path=run_work_dir
            / _Files.PARTICLES_VS_ZONES_ASCII,
        )
        # b) Parse zones in voids raw file in the work directory
        _postprocessing.parse_zones_in_void_output(
            executable_path=_Paths.ZOBOV
            / _ExecutableNames.ZONES_IN_VOIDS_BIN,
            input_file_path=run_work_dir
            / _Files.ZONES_VS_VOID_RAW,
            output_file_path=run_work_dir
            / _Files.ZONES_VS_VOID_ASCII,
        )
        # Process 2:
        # a) Get Tuple of (VoidProperties, particles)
        zobov_vp_and_part = (
            _postprocessing.process_and_extract_void_properties_and_particles(
                tinz_output_file_path=run_work_dir
            / _Files.PARTICLES_VS_ZONES_ASCII,
                zinv_output_file_path=run_work_dir
            / _Files.ZONES_VS_VOID_ASCII,
                jozov_text_file_output_path=run_work_dir
            / _Files.OUTPUT_JOZOV_VOIDS_DAT
            )
            )

        # b) divide the output
        particle_by_voids, zobov_voids = [], []
        for void_properties, particle_in_void in zobov_vp_and_part:
            particle_by_voids.append(particle_in_void)
            zobov_voids.append(void_properties)

        # c) Create extra
        extra = {
            "zobov_path": self._zobov_path,
            "zobov_voids": tuple(zobov_voids),
            "files_dir": run_work_dir
        }

        return tuple(particle_by_voids), extra
