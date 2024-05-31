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
path_src = path.Path(os.path.abspath("./voidfindertk/zobov/src"))

params = {
    "path_src": path_src/"src",
    "path_input_file": path_src/"src"/"tracers_zobov.raw",
    "buffer_size":0.08,
    "box_size":500,
    "number_of_divisions": 2,
    "executable_name":"output_vozinit",
    "output_name_particles_in_zones":"part_vs_zone",
    "output_name_zones_in_void":"zones_vs_voids",
    "output_name_text_file":"ouput_txt",
    "density_threshold": 0
    }

wrap.run_zobov_void_finder(**params)
