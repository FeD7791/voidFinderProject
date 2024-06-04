"""
Wrapper functions for the ZOBOV void finder.

This module provides wrapper executables to run the ZOBOV void finder algorithm
and its associated preprocessing steps.

"""

import ctypes
import os
import shutil

import numpy as np

import sh

from ..utils import chdir


def _move_inputs(src, dst_dir):
    """Move input files to the destination directory if necessary.

    Parameters:
    -----------
    src : str
        Path to the source file.
    dst_dir : str
        Path to the destination directory.

    Returns:
    --------
    str
        Path to the moved or original file.
    """
    src_dir = os.path.dirname(src)
    if src_dir == dst_dir:
        return src
    return shutil.copy2(src, dst_dir)


def run_vozinit(
    *,
    vozinit_dir_path,
    input_file_path,
    buffer_size,
    box_size,
    number_of_divisions,
    executable_name,
    work_dir_path,
):
    """
    Run the vozinit command of the ZOBOV void finder.

    Parameters:
    -----------
    vozinit_dir_path : str
        Path to the directory containing the vozinit executable.
    input_file_path : str
        Path to the tracers input file.
    buffer_size : float
        Sets the size, in units such that the box size of the data cube is 1.
    box_size : float
        The range of positions of particles in each dimension.
    number_of_divisions : int
        Number of divisions (default 2)
        The no. of partitions in each dimension; must be at least 2.
    executable_name : str
        Suffix label for the output files
    work_dir_path : str
        Path to the working directory.

    Returns:
    --------
    str
        Output of the vozinit command.
    """

    vozinit = sh.Command("vozinit", search_paths=[vozinit_dir_path])

    params = (
        input_file_path,
        buffer_size,
        box_size,
        number_of_divisions,
        executable_name,
    )

    args = [str(param) for param in params]

    with chdir(work_dir_path):
        output = vozinit(*args)

    return output


def run_voz_step(
    *,
    preprocess_dir_path,
    executable_name,
    work_dir_path,
    voz_executables_path,
):
    """
    Run the preprocessing step of the ZOBOV void finder.

    Parameters:
    -----------
    preprocess_dir_path : str
        Path to the directory containing the preprocessing executable.
    executable_name : str
        Suffix of the executable file.
    work_dir_path : str
        Path to the working directory.
    voz_executables_path : str
        Path of the voz - 1b1 / tie (voz1b1/voztie) exectutables

    Returns:
    --------
    str
        Output of the preprocessing command.
    """
    # Moving necesary files tho run the src-executable
    _move_inputs(voz_executables_path / "voz1b1", work_dir_path)
    _move_inputs(voz_executables_path / "voztie", work_dir_path)

    full_executable_name = preprocess_dir_path / f"scr{executable_name}"
    preprocess = sh.Command(full_executable_name)

    with chdir(work_dir_path):
        output = preprocess()

    return output


def run_voz1b1(
    *,
    input_file_path,
    buffer_size,
    box_size,
    executable_name,
    number_of_divisions,
    binary_division,
    voz1b1_dir_path,
    work_dir_path,
):
    """
    Run the VOZ1B1 executable with specified parameters.

    Parameters:
    ----------
    input_file_path : str
        Path to the input file.
    buffer_size : int
        Sets the size, in units such that the box size of the data cube is 1.
    box_size : int
        Max value of the tracers in each x,y,z coordinates.
    executable_name : str
        Suffix of the created executable.
    number_of_divisions : int
        The no. of partitions in each dimension;
        Must be at least 2 (giving 8 sub-boxes)
    binary_division : tuple
        Tuple of 3 integers representing binary divisions.
    voz1b1_dir_path : str
        Path to the voz1b1 executable directory.
    work_dir_path : str
        Path to the working directory.

    Returns:
    -------
    str
        Output of the command.
    """

    voz1b1 = sh.Command(voz1b1_dir_path / "voz1b1")

    params = [
        input_file_path,
        buffer_size,
        box_size,
        executable_name,
        number_of_divisions,
        binary_division[0],
        binary_division[1],
        binary_division[2],
    ]

    args = [str(param) for param in params]
    with chdir(work_dir_path):
        output = voz1b1(*args)
    return output


def run_voztie(
    *,
    number_of_divisions,
    executable_name,
    voztie_dir_path,
    work_dir_path,
):
    """
    Run the VOZTIE executable with specified parameters.

    Parameters:
    ----------
    number_of_divisions : int
        Number of divisions, min: 2.
    executable_name : str
        Suffix of the executable.
    voztie_dir_path : str
        Path to the VOZTIE directory.
    work_dir_path : str
        Path to the working directory.

    Returns:
    -------
    str
        Output of the sh voztie Command.
    """
    voztie = sh.Command(voztie_dir_path / "voztie")

    params = [
        number_of_divisions,
        executable_name,
    ]
    args = [str(param) for param in params]
    with chdir(work_dir_path):
        output = voztie(*args)
    return output


def run_jozov(
    *,
    jozov_dir_path,
    executable_name,
    output_name_particles_in_zones,
    output_name_zones_in_void,
    output_name_text_file,
    density_threshold,
    work_dir_path,
):
    """
    Run the jozov command of the ZOBOV void finder.

    Parameters:
    -----------
    jozov_dir_path : str
        Path to the directory containing the jozov executable.
    executable_name : str
        Corpus Name of the executable files with suffix adj/vol.
        These are the adjacency and volume input files of the run.
    output_name_particles_in_zones : str
        Name of the output file containing particles vs zones.
    output_name_zones_in_void : str
        Name of the output file containing zones vs voids.
    output_name_text_file : str
        Name of the main parameters output text file.
    density_threshold : float
        Density threshold (0 for no threshold).
    work_dir_path : str
        Path to the working directory.

    Returns:
    --------
    str
        Output of the jozov command.
    """
    jozov = sh.Command(jozov_dir_path / "jozov")

    args = (
        f"adj{executable_name}.dat",
        f"vol{executable_name}.dat",
        f"{output_name_particles_in_zones}.dat",
        f"{output_name_zones_in_void}.dat",
        f"{output_name_text_file}.dat",
        f"{density_threshold}",
    )

    with chdir(work_dir_path):
        output = jozov(*args)

    return output


def write_input(
    box, path_executable, path_raw_file_output, path_txt_file_output
):
    """
    Create input binary files for the Zobov finder using an input
    Box of tracers.

    Parameters:
    ----------
    box : object
        Input box object containing x, y, z, vx, vy, vz, and m.
    path_executable : str
        Path to the executable file.
    path_raw_file_output : str
        Path to the raw file output directory.
    path_txt_file_output : str
        Path to the text file output directory.

    Returns:
    -------
    None
    """

    # Create input binary files for Zobov finder

    # Declare library path
    clibrary = ctypes.CDLL(str(path_executable), mode=ctypes.RTLD_GLOBAL)

    # Create Input Pointers for x,y,z,vx,vz,vy,m
    arr_pointer = 7 * [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=["CONTIGUOUS"])
    ]

    # Declare Input Pointers type
    clibrary.c_binary_writter.argtypes = arr_pointer + [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    # Fill Input
    clibrary.c_binary_writter(
        box.x,
        box.y,
        box.z,
        box.vx,
        box.vy,
        box.vz,
        box.m,
        len(box),
        os.path.join(path_raw_file_output, "tracers_zobov.raw").encode(
            "utf-8"
        ),
        os.path.join(path_txt_file_output, "tracers_zobov.txt").encode(
            "utf-8"
        ),
    )
