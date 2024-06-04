import os
import pathlib as path
import tempfile

import numpy as np

import sh

from voidfindertk.io import read_table
from voidfindertk.zobov import _wrapper as wrap

path_dataset = "./datasets/halos_fede.dat"
path_loader_so = "./voidfindertk/zobov/src/zobov_loader.so"


box = read_table(path_dataset)


path_src = path.Path(os.path.abspath("./voidfindertk/zobov/src"))


with tempfile.TemporaryDirectory(suffix="_zovob", dir=path_src) as work_dir:
    ##
    work_dir = path.Path(work_dir)
    ##
    wrap.write_input(
        box=box,
        path_executable=path_loader_so,
        path_raw_file_output=work_dir,
        path_txt_file_output=work_dir,
    )

    voz1b1_params = {
        "voz1b1_dir_path": path_src / "src",
        "input_file_path": work_dir / "tracers_zobov.raw",
        "buffer_size": 0.08,
        "box_size": 500,
        "executable_name": "output_vozinit",
        "number_of_divisions": 2,
        "binary_division": [0, 0, 0],
        "work_dir_path": work_dir,
    }

    for i in np.arange(8):
        voz1b1_params["binary_division"] = [*format(i, "03b")]
        output_voz1b1 = wrap.run_voz1b1(**voz1b1_params)

    output_voztie = wrap.run_voztie(
        work_dir_path=work_dir,
        voztie_dir_path=path_src / "src",
        number_of_divisions=2,
        executable_name="output_vozinit",
    )

    output_jozov = wrap.run_jozov(
        jozov_dir_path=path_src / "src",
        executable_name="output_vozinit",
        output_name_particles_in_zones="part_vs_zone",
        output_name_zones_in_void="zones_vs_voids",
        output_name_text_file="output_txt",
        density_threshold=0,
        work_dir_path=work_dir,
    )
    print(sh.ls(work_dir))
