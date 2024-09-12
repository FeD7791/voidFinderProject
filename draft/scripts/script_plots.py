from voidfindertk import read_table
import pathlib
from voidfindertk.core import vsf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dataset_path = "/home/jorgefederico/updates/vftk_actual002/voidFinderProject/datasets/halos_ascii_1000_1024_npmin_10_z0.51.dat"
workdir_path = pathlib.Path("/home/jorgefederico/updates/vftk_actual002/voidFinderProject/run_popcorn")

dbox = read_table(dataset_path,names=["m","x", "y", "z", "vx", "vy", "vz"])

dir = "tmplxhn4wcw2024-08-22T01:08:51.465883+00:00"
df = pd.read_csv(
    str(workdir_path/f"{dir}"/"sphfile.dat"),
    names=["id","rad","x","y","z","delta"],
    delim_whitespace=True
    )
rad = np.array(df["rad"])


voidsizefun = vsf.void_size_function(
    radius=rad,
    box=dbox.box,
    delta=-0.9,
    scale_2_num_samples=1
)

param1 = "hola"
param2 = "mundo"
plt.figure(figsize=(8, 8))
plt.plot(mids, normalized_counts, marker='o', linestyle='-', label='Void Size Function')
plt.xlabel(r'$log_{10}(R)$ '+ f' {dbox.box.x[0].unit}')
plt.ylabel(r'$\frac{1}{V} \frac{dN_v}{dlnR_v}$', fontsize=20)
# plt.xlim((mids[np.where(mids>0)[0][0]],mids[-1]))
plt.yscale('log')
plt.title('Void Size Function')
plt.grid(True)
plt.legend()
plt.text(min(mids),min(normalized_counts), f"{param1}\n {param2}", fontsize = 22, 
         bbox = dict(facecolor = 'cyan', alpha = 0.5))
plt.savefig('void_size_function43.jpg')
# plt.show()



#./clean_duplicates config=/home/jorgefederico/updates/vftk_actual002/voidFinderProject/run_popcorn/tmplxhn4wcw2024-08-22T01:08:51.465883+00:00/vars.conf