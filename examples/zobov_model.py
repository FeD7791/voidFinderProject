from voidfindertk import zobov
import os
import pathlib
import tempfile

from voidfindertk.io import read_table

# Path to this file
# This path is used to ensure that example works anywhere
this_file_path = pathlib.Path(os.path.abspath(__file__))

# Path to the dataset file
this_file_parent = this_file_path.parents[1]
# path_dataset = this_file_parent/'datasets'/'halos_fede.dat'
path_dataset = "./voidFinderProject/datasets/halos_fede.dat"

# Generate Directory to place the runs
temp_dir = tempfile.mkdtemp(prefix="run", dir=str(this_file_parent))
# print(str(pathlib.Path(temp_dir)))

# Generate the model
zobov_0 = zobov.ZobovVF(
    # Perform zobov and place its files in a directory inside temp_dir
    workdir=str(str(pathlib.Path(temp_dir)))
)

# Get dataset box object:
tracers = read_table.read_table(path_or_buffer=path_dataset)

# Perform Void Finder method on dataset
voids = zobov_0.find(tracers)

# #workdir
# work_dir = pathlib.Path("/home/jorgefederico/updates/test_vftk/voidFinderProject/runp3o9zuay/tmpr4ni14vn2024-07-22T00:14:26.827546+00:00")
# #zobov_path
# zobov_path = pathlib.Path("/home/jorgefederico/updates/test_vftk/testenv/lib/python3.11/site-packages/voidfindertk/zobov/src")


# def mockfun():
#     zobov_properties = zobov.process_and_extract_void_properties_and_particles(
#     zinv_executable_path=zobov_path/"zones_in_void.so",
#     tinz_executable_path=zobov_path/"tracers_in_zones.so",

#     tinz_input_file_path=work_dir/"part_vs_zone.dat",
#     zinv_input_file_path=work_dir/"zones_vs_voids.dat",

#     tinz_output_file_path=work_dir/"part_vs_zone.txt",
#     zinv_output_file_path=work_dir/"zones_vs_voids.txt",

#     jozov_text_file_output_path=work_dir/"output_txt.dat"
#     )

#     particle_by_voids, zobov_voids = [], []
#     for void_properties, particle_in_void in zobov_properties:
#         particle_by_voids.append(particle_in_void)
#         zobov_voids.append(void_properties)

#     extra = {
#                 "zobov_voids": tuple(zobov_voids),
#             }
#     return tuple(particle_by_voids),extra

# data1,data2 = mockfun()

# voids = Voids(
#     method="voidcito",
#     tracers=databox.box,
#     voids=voids_tuple,
#     extra=extra,
# )
