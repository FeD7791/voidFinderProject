from voidfindertk import zobov, io
import pathlib
import os

dataset_path = pathlib.Path("./datasets")
workdir_path = pathlib.Path(os.path.abspath("./runz"))

box = io.read_table(
    dataset_path/"halos_fede.dat",
    names=["x", "y", "z", "vx", "vy", "vz","m"])

model = zobov.ZobovVF(
    density_threshold=0.1, box_size=500, workdir=workdir_path
)
# run_work_dir = pathlib.Path('/home/jorgefederico/updates/vftk_1109/voidFinderProject/runz/tmpqt51ba9_2024-09-23T19:33:00.395859+00:00')
# model_find_parameters = {"run_work_dir": run_work_dir, "box": box}
# void = model.build_voids(model_find_parameters=model_find_parameters)
void = model.find(box)
#joblib.dump(output, "/home/jorgefederico/updates/vftk_1109/voidFinderProject/tests/mock_data/run_jozov.jpkl")



grid = gsp.GriSPy(cloud_with_voids, copy_data=False, N_cells=64)
periodic = {0:(0,1000),1:(0,1000),2:(0,1000)}
grid.set_periodicity(periodic, inplace=True)
lev_dist, lev_ind = grid.nearest_neighbors(centers[14:15],100)

ro = 3*np.arange(1,101)/(4*np.pi*lev_dist[0]**3)
dist = (3*np.arange(1,len(d0[14])+1)/(4*np.pi*d0[14]))**(1/3)

dist[np.where(d0[14]<rho_threshold)[0][-1]]

import matplotlib.pyplot as plt
plt.scatter(np.arange(1,len(d0[14])+1),d0[14])