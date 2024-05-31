from voidfindertk.zobov import _wrapper as wrap
from voidfindertk.zobov._zobov import write_zobov_input
from voidfindertk.io import read_table
import pathlib as path
import os


path_dataset = "./datasets/halos_fede.dat"
path_loader_so = "./voidfindertk/zobov/src/zobov_loader.so"
path_zobov_src_src = "./voidfindertk/zobov/src/src"

# box = read_table(path_dataset)

# write_zobov_input(
#     box=box,
#     path_executable=path_loader_so,
#     path_raw_file_output=path_zobov_src_src,
#     path_txt_file_output=path_zobov_src_src
#     )

##############################################################
work_dir = None

import tempfile

path_src = path.Path(os.path.abspath("./voidfindertk/zobov/src"))


with tempfile.TemporaryDirectory(suffix="_zovob") as work_dir:

    output_vozinit = wrap.run_vozinit(
        vozinit_dir_path=path_src / "src",
        input_file_path=path_src / "src" / "tracers_zobov.raw",
        buffer_size=0.08,
        box_size=500,
        number_of_divisions=2,
        executable_name="output_vozinit",
        output_name_particles_in_zones="part_vs_zone",
        output_name_zones_in_void="zones_vs_voids",
        output_name_text_file="ouput_txt",
        density_threshold=0,
        work_dir_path=work_dir,
    )

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
        output_name_text_file="ouput_txt",
        density_threshold=0,
        work_dir_path=work_dir,
    )
