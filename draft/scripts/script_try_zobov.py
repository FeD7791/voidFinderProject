from voidfindertk import io
from voidfindertk import zobov
from voidfindertk.voids import Voids
import time
import pathlib

# Run full finder

path = "/home/jorgefederico/updates/vftk_actual002/voidFinderProject/datasets/halos_ascii_1000_1024_npmin_10_z0.51.dat"
# Get Box
dbox = io.xyz_read_table(path, usecols=[1, 2, 3])
box = dbox.box

workdir = '/home/jorgefederico/updates/vftk_1109/voidFinderProject/runz'

# model
model = zobov.ZobovVF(box_size=1000, workdir=workdir)

# start_time = time.time()
# voids = model.find(databox=dbox)
# print("--- %s seconds ---" % (time.time() - start_time))

## Just build voids
model_find_parameters = {
    "run_work_dir": pathlib.Path(workdir) / "old",
    "box": box,
}

start_time = time.time()
tinv, centers, extra = model.build_voids(
    model_find_parameters=model_find_parameters
)
##
voids = Voids(
    method="zobov",
    box=box,
    tracers_in_voids=tinv,
    centers=centers,
    extra=extra,
)
errors, rad, tracers, density = voids.all_effective_radius()
print("--- %s seconds ---" % (time.time() - start_time))
