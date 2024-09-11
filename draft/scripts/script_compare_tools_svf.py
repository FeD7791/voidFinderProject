from voidfindertk.svf_popcorn import _svf_popcorn
from voidfindertk import io,tools
import os
import pathlib
import pandas as pd
import numpy as np


dataset_path = pathlib.Path("./datasets")
workdir_path = pathlib.Path(os.path.abspath("/home/jorgefederico/updates/vftk_actual002/voidFinderProject/run_popcorn"))
# run_workdir_path = workdir_path / "tmpovse6ej82024-08-21T16:12:49.132243+00:00" #-0.8
# run_workdir_path = workdir_path / "tmpw990bohs2024-08-09T18:06:26.180577+00:00" #-0.9
run_workdir_path = workdir_path / "tmplxhn4wcw2024-08-22T01:08:51.465883+00:00" #-0.7

dbox = io.read_table(dataset_path/"halos_ascii_1000_1024_npmin_10_z0.51.dat",names=["m","x", "y", "z", "vx", "vy", "vz"])
#Params
delta = -0.9
n_neighbors = 1000

model = _svf_popcorn.PopCornVF(
    workdir=workdir_path,
    boxsize=1000,
    densth=delta,
    minradius = 5,
    maxradius = 100
    )

##Just run build_voids
model_find_parameters = {
    "run_work_dir": run_workdir_path,
    "box": dbox.box}
tinv,centers,extra = model.build_voids(model_find_parameters)

df = pd.read_csv(run_workdir_path / "sphfile.dat",delim_whitespace=True, names=['ID','Rad','x','y','z','delta'])
centers = np.array(df[['x','y','z']])
rad,tracers,density = tools.calculate_r_eff(centers=centers,box=dbox.box,delta=delta,n_neighbors=n_neighbors)
rad2 = np.array(rad)






n_nat = np.arange(1,n_neighbors+1)
crit_density = (1+delta)*(len(dbox.box)/(dbox.box.size()**3))
average_density = len(dbox.box)/(dbox.box.size()**3)

#rads type -2
len(np.where(rad2 == -2)[0])
#rads type 0
len(np.where(rad2 == 0)[0])
#rads type -1
len(np.where(rad2 == -1)[0])
#Void that are not type -1 or 0
old_index = np.arange(0,len(rad))
#Indices de voids que no son ni tipo -1 ni tipo -2
index = np.where((rad2 != 0) & (rad2 !=-1) & (rad2 !=-1))[0]

abs_dist = abs(np.array(df['Rad'][index]) - rad2[index])
# Percentage
100*len(np.where(abs_dist<1)[0])/len(abs_dist)

# fig = plt.figure(figsize=(7,7))
# ax = fig.add_subplot(1,1,1)
# ax.scatter(n_nat,density[23944], s=0.5)
# ax.axhline(y=average_density, label=f"{average_density:.4f}",c='green',ls='--')
# ax.axhline(y=crit_density, label=f"{crit_density:.4f}", c='red',ls='--')
# ax.set_xlabel("N tracers")
# ax.set_ylabel("Density")
# ax.legend()
# plt.savefig("rad_dens.jpg")



fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(1,1,1)
ax.scatter(n_nat,c[2], s=0.5)
ax.axhline(y=average_density, label=f"{average_density:.4f}",c='green',ls='--')
ax.axhline(y=crit_density, label=f"{crit_density:.4f}", c='red',ls='--')
ax.set_xlabel("N tracers")
ax.set_ylabel("Density")
ax.legend()
plt.savefig("rad_dens.jpg")














