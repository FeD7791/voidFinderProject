# This script is buil to find the reff of zobov voids and analyze them

from voidfindertk import io, tools
import os
import pathlib
import pandas as pd
import numpy as np

dataset_path = pathlib.Path("./datasets")
workdir_path = pathlib.Path(
    os.path.abspath(
        "/home/jorgefederico/updates/vftk_actual002/voidFinderProject/runv"
    )
)
run_workdir_path = (
    workdir_path / "tmpaizyo2s62024-09-07T21:10:14.262424+00:00"
)  # 0.8
run_workdir_path = (
    workdir_path / "tmpaizyo2s62024-09-07T21:10:14.262424+00:00"
)  # 0.9
# run_workdir_path = workdir_path / "tmpzcxk8adm2024-09-07T21:58:35.096221+00:00"#0.7
dbox = io.read_table(
    dataset_path / "halos_ascii_1000_1024_npmin_10_z0.51.dat",
    names=["m", "x", "y", "z", "vx", "vy", "vz"],
)
box = dbox.box

df = pd.read_csv(
    run_workdir_path / "output_txt.dat", delim_whitespace=True, header=1
)
center_index = np.array(df["CoreParticle"])
xyz = np.column_stack((box.x.value, box.y.value, box.z.value))
centers = xyz[center_index]
n_neighbors = 400
delta = -0.9

rad, tracers, density = tools.calculate_r_eff(
    centers=centers, box=dbox.box, delta=delta, n_neighbors=n_neighbors
)
rad2 = np.array(rad)

n_nat = np.arange(1, n_neighbors + 1)
crit_density = (1 + delta) * (len(dbox.box) / (dbox.box.size() ** 3))
average_density = len(dbox.box) / (dbox.box.size() ** 3)

# rads type -2
len(np.where(rad2 == -2)[0])
# rads type 0
len(np.where(rad2 == 0)[0])
# rads type -1
len(np.where(rad2 == -1)[0])
# rads normal
len(np.where((rad2 != 0) & (rad2 != -1) & (rad2 != -2))[0])
# Void that are not type -1 or 0
old_index = np.arange(0, len(rad))
index = np.where((rad2 != 0) & (rad2 != -1) & (rad2 != -1))[0]

max(rad2[np.where((rad2 != 0) & (rad2 != -1) & (rad2 != -2))[0]])
