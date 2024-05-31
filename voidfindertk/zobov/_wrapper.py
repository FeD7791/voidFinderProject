import os
import shutil

import sh

from ..utils import chdir


def _move_inputs(src, dst_dir):
    src_dir = os.path.dirname(src)
    if src_dir == dst_dir:
        return src
    return shutil.copy2(src, dst_dir)


# =============================================================================
# RUN ZOBOV
# =============================================================================


def run_vozinit(
    *,
    vozinit_dir_path,
    input_file_path,
    buffer_size,
    box_size,
    number_of_divisions,
    executable_name,
    output_name_particles_in_zones,
    output_name_zones_in_void,
    output_name_text_file,
    density_threshold,
    work_dir_path,
):

    input_file_path = _move_inputs(input_file_path)

    # Runing Zobov
    vozinit = sh.Command(
        "vozinit", search_paths=[vozinit_dir_path]
    )  # path_src --> pathlib.Path

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
    full_executable_name = preprocess_dir_path / f"scr{executable_name}"
    preprocess = sh.Command(full_executable_name)

    with chdir(work_dir_path):
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
    jozov = sh.Command(jozov_dir_path / "jozov")

    args = (
        f"adj{executable_name}.dat",
        f"vol{executable_name}",
        f"{output_name_particles_in_zones}.dat",
        f"{output_name_zones_in_void}.dat",
        f"{output_name_text_file}.dat",
        f"{density_threshold}",
    )

    with chdir(work_dir_path):
        output = jozov(*args)

    return output
