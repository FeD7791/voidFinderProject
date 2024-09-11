from voidfindertk import io
from voidfindertk import dive
import time
import pathlib

# Run full finder

path = "/home/jorgefederico/updates/vftk_actual002/voidFinderProject/datasets/halos_ascii_1000_1024_npmin_10_z0.51.dat"
# path = "/home/jorgefederico/updates/vftk_actual002/voidFinderProject/datasets/halo_catalogue.txt"
# Get Box
dbox = io.xyz_read_table(path, usecols=[1,2,3])
box = dbox.box

workdir = "/home/jorgefederico/updates/vftk_actual002/voidFinderProject/runv"

# model
model = dive.DiveVF(box_size=1000, workdir=workdir, delta_r=[0.,50.], overlap_criterion=True, ratio=0.1)
# start_time = time.time()
# voids = model.find(databox=dbox)
# print("--- %s seconds ---" % (time.time() - start_time))

## Just build voids
model_find_parameters = {
    'run_work_dir':pathlib.Path(workdir) / 'tmpluh7pb6n2024-08-09T19:14:07.847469+00:00',
    'box':box
                         }

start_time = time.time()
v1,v2,v3 = model.build_voids(model_find_parameters=model_find_parameters)
print("--- %s seconds ---" % (time.time() - start_time))