from voidfindertk import io, popcorn
import os
import pathlib
import time

dataset_path = pathlib.Path("./datasets")
workdir_path = pathlib.Path(os.path.abspath("./run_real_popcorn"))

dbox = io.read_table(dataset_path/"halos_ascii_1000_1024_npmin_10_z0.51.dat",names=["m","x", "y", "z", "vx", "vy", "vz"])

model = popcorn.PopCorn(
    workdir=workdir_path,
    boxsize=1000,
    densth=-0.9,
    minradius = 20,
    maxradius = 100,
    shot_noise_threshold = 25
)
void = model.find(databox=dbox)



# read_and_modify_config(
#     config_file_path="./vars.conf",
#     section="INPUT_PARAMS",
#     parameter="MINRADIUS",
#     new_value=str(57)
#     )