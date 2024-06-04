"""
Wrapper functions for the ZOBOV void finder.

This module provides wrapper executables to run the ZOBOV void finder algorithm
and its associated preprocessing steps.

"""

import os
import shutil

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
        Path to the input file.
    buffer_size : float
        Buffer size parameter for vozinit.
    box_size : float
        Box size parameter for vozinit.
    number_of_divisions : int
        Number of divisions parameter for vozinit.
    executable_name : str
        Name of the executable file.
    output_name_particles_in_zones : str
        Name of the output file for particles in zones.
    output_name_zones_in_void : str
        Name of the output file for zones in voids.
    output_name_text_file : str
        Name of the output text file.
    density_threshold : float
        Density threshold parameter for vozinit.
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


def run_preprocess(*, preprocess_dir_path, executable_name, work_dir_path):
    """Run the preprocessing step of the ZOBOV void finder.

    Parameters:
    -----------
    preprocess_dir_path : str
        Path to the directory containing the preprocessing executable.
    executable_name : str
        Name of the executable file.
    work_dir_path : str
        Path to the working directory.

    Returns:
    --------
    str
        Output of the preprocessing command.
    """
    full_executable_name = preprocess_dir_path / f"scr{executable_name}"
    preprocess = sh.Command(full_executable_name)

    with chdir(work_dir_path):
        print(os.getcwd())
        output = preprocess()

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
    """Run the jozov command of the ZOBOV void finder.

    Parameters:
    -----------
    jozov_dir_path : str
        Path to the directory containing the jozov executable.
    executable_name : str
        Name of the executable file.
    output_name_particles_in_zones : str
        Name of the output file for particles in zones.
    output_name_zones_in_void : str
        Name of the output file for zones in voids.
    output_name_text_file : str
        Name of the output text file.
    density_threshold : float
        Density threshold parameter for jozov.
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
