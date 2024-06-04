import os
import pathlib as path
import tempfile

import sh

from voidfindertk.io import read_table
from voidfindertk.zobov import _wrapper as wrap
from voidfindertk.zobov._zobov import write_zobov_input

path_dataset = "./datasets/halos_fede.dat"
path_loader_so = "./voidfindertk/zobov/src/zobov_loader.so"


box = read_table(path_dataset)


path_src = path.Path(os.path.abspath("./voidfindertk/zobov/src"))


with tempfile.TemporaryDirectory(suffix="_zovob", dir=path_src) as work_dir:
    ##
    work_dir = path.Path(work_dir)
    ##
    write_zobov_input(
        box=box,
        path_executable=path_loader_so,
        path_raw_file_output=work_dir,
        path_txt_file_output=work_dir,
    )

    output_vozinit = wrap.run_vozinit(
        vozinit_dir_path=path_src / "src",
        input_file_path=work_dir / "tracers_zobov.raw",
        buffer_size=0.08,
        box_size=500,
        number_of_divisions=2,
        executable_name="output_vozinit",
        work_dir_path=work_dir,
    )

    wrap._move_inputs(path_src / "src" / "voz1b1", work_dir)
    wrap._move_inputs(path_src / "src" / "voztie", work_dir)

    output_preprocess = wrap.run_preprocess(
        preprocess_dir_path=work_dir,
        executable_name="output_vozinit",
        work_dir_path=work_dir,
    )

    wrap.run_jozov(
        jozov_dir_path=path_src / "src",
        executable_name="output_vozinit",
        output_name_particles_in_zones="part_vs_zone",
        output_name_zones_in_void="zones_vs_voids",
        output_name_text_file="output_txt",
        density_threshold=0,
        work_dir_path=work_dir,
    )
