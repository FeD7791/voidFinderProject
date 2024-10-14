from voidfindertk.svf_popcorn import PopCornVF
from voidfindertk import read_table
from voidfindertk.core import Voids
import os
import pathlib
import time

dataset_path = pathlib.Path("./datasets")
workdir_path = pathlib.Path(os.path.abspath("./run_popcorn"))

dbox = read_table(
    dataset_path / "halos_ascii_1000_1024_npmin_10_z0.51.dat",
    names=["m", "x", "y", "z", "vx", "vy", "vz"],
)

start_time = time.time()
model = PopCornVF(
    workdir=workdir_path,
    boxsize=1000,
    densth=-0.7,
    minradius=20,
    maxradius=100,
)
# void = model.find(databox=dbox)

##Just run build_voids
model_find_parameters = {
    "run_work_dir": workdir_path
    / "tmp4ut8qrcf2024-08-22T01:49:27.150943+00:00",
    "box": dbox.box,
}
tinv, centers, extra = model.build_voids(model_find_parameters)
##
voids = Voids(
    method="popcorn",
    box=dbox.box,
    tracers_in_voids=tinv,
    centers=centers,
    extra=extra,
)
errors, rad, tracers, density = voids.effective_radius()
objeto = voids.effective_radius()
print("--- %s seconds ---" % (time.time() - start_time))
