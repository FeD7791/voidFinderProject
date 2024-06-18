import datetime as dt
import os
import pathlib
import shutil
import tempfile

import numpy as np

from ..models import ModelABC
from . import _wrapper as _wrap


class _Paths:
    """
    Class that holds paths of reference to the current file and
    ZOBOV's src directory
    """
    CURRENT = pathlib.Path(
        os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    )
    ZOBOV = CURRENT / "src"  # Path to the src folder of Zobov


class _Files:
    """
    
    """
    TRACERS_RAW = "tracers_zobov.raw"
    TRACERS_TXT = "tracers_zobov.txt"


class ZobovVF(ModelABC):
    """
    ZobovVF class for running ZOBOV Void Finder.

    This class provides methods to preprocess data and execute the ZOBOV
    Void Finder algorithm on a given data box.

    Parameters
    ----------
    buffer_size : float, optional
        Buffer size for ZOBOV (default is 0.08).
    box_size : int, optional
        Size of the box containing the data (default is 500).
    number_of_divisions : int, optional
        Number of divisions in each dimension of the box (default is 2).
    density_threshold : int, optional
        Density threshold for void identification (default is 0).
    zobov_path : str or None, optional
        Path to ZOBOV executable (default is None, uses internal path).
    workdir : str or None, optional
        Temporary working directory path (default is None, creates a new temp directory).
    workdir_clean : bool, optional
        Whether to clean up the working directory on deletion (default is False).
    dtype : numpy.dtype, optional
        Data type used for computations (default is np.float32).

    Methods
    -------
    preprocess(databox)
        Placeholder method for data preprocessing.
    model_find(databox)
        Executes the ZOBOV Void Finder algorithm on the provided DataBox object.

    Notes
    -----
    The ZOBOV Void Finder is executed in several steps including VOZINIT, VOZSTEP, and JOZOV.
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
        self._ensity_threshold = density_threshold

        self._zobov_path = pathlib.Path(
            _Paths.ZOBOV if zobov_path is None else zobov_path
        )
        # Create a workdir path to run ZOBOV
        self._workdir = pathlib.Path(
            tempfile.mkdtemp(prefix=f"vftk_{type(self).__name__}_")
            if workdir is None
            else workdir
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
        return self._ensity_threshold

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
        This method will create a temporal directory inside the working directory 
        of the ZobovVF class workdir.

        Returns
        -------
            run_work_dir(pathlib.Path): path of the work directoty 
        """
        timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
        run_work_dir = pathlib.Path(
            tempfile.mkdtemp(suffix=timestamp, dir=self.workdir)
        )
        return run_work_dir
    
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
        Execute the ZOBOV Void Finder algorithm on the provided DataBox object.

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
            path_executable=self._zobov_path / "zobov_loader.so",
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
            executable_name="output_vozinit",
            work_dir_path=run_work_dir,
        )

        # VOZSTEP =============================================================
        # This step is mandatory if VOZINIT was run before

        _wrap.run_voz_step(
            preprocess_dir_path=run_work_dir,
            executable_name="output_vozinit",
            work_dir_path=run_work_dir,
            voz_executables_path=_Paths.ZOBOV
            / "src",  # this is the path where voz1b1 and voztie exe are
        )

        # JOZOV =============================================================
        _wrap.run_jozov(
            jozov_dir_path=_Paths.ZOBOV / "src",
            executable_name="output_vozinit",
            output_name_particles_in_zones="part_vs_zone",
            output_name_zones_in_void="zones_vs_voids",
            output_name_text_file="output_txt",
            density_threshold=0,
            work_dir_path=run_work_dir,
        )

    def __del__(self):
        """
        Destructor that cleans up the temporary working directory 
        if workdir_clean is True.
        """
        if self._workdir_clean:
            shutil.rmtree(self._workdir)
