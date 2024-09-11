import grispy as gsp
import numpy as np
import matplotlib.pyplot as plt
import pathlib
import pandas as pd
from voidfindertk import io
# find tracers that have an amount of density

dataset = pathlib.Path("./datasets/") / 'halos_ascii_1000_1024_npmin_10_z0.51.dat'

#zobov_centers = pathlib.Path("./runz")
#dbox = io.read_table(data)


dbox = io.read_table(dataset,names=['m','x','y','z','vx','vy','vz'])
box = dbox.box
xyz = np.column_stack((box.x.value,box.y.value,box.z.value))

grid = gsp.GriSPy(xyz)
###zobov
# df = pd.read_csv(zobov_centers/"tmpe7gbr1792024-08-04T18:23:45.084251+00:00"/"output_txt.dat",sep='\s+',header=1)
# centers = xyz[df['CoreParticle']]
run_pocorn = pathlib.Path('/home/jorgefederico/updates/vftk_actual002/voidFinderProject/run_popcorn')
sphfile = run_pocorn / 'tmplxhn4wcw2024-08-22T01:08:51.465883+00:00' / 'sphfile.dat'

df = pd.read_csv(sphfile,delim_whitespace=True, names=['ID','Rad','x','y','z','delta'])
centers = np.array(df[['x','y','z']])

# Esto me da para cada centro, la distancia de los n puntos de la grilla mas cercanos a este y su indice
dist,nn=grid.nearest_neighbors(centres=centers,n=100)

tracers_new = []
rad_new = []
density_values = []
n_nat = np.arange(len(dist[0]))
crit_density = 0.3*(len(box)/(box.size()**3))
for n,d in enumerate(dist):
    # Find density values for n_nat particles at radius d
    density_n_nat_d = (3*n_nat[1:])/(4*np.pi*d[1:]**3)  ######Starts with index 0
    # Find all density values that are less than crit_density
    dens_values = np.where(density_n_nat_d <crit_density)[0] ##### Los indices quehay aqui empiezan con 0
    density_values.append(density_n_nat_d)
    if(len(dens_values)==0):
        pass
    else:
        # From the values that fulfill the latter condition find the index of the
        # value with max radii
        dist_max_index = np.where(d[dens_values+1]==max(d[dens_values+1]))[0][0]

        # Give the all the tracers whose distance is lesser than the found radii
        tracers_new.append(nn[n][:dist_max_index])
        # Final radii is half distance between distk_max_index and dist_max_index +1
        try:
            rad_new.append((d[dist_max_index+1]+d[dist_max_index])/2)
        except IndexError:
            print(f"problem with value{n}")
        # rad_new.append(d[dist_max_index])


#plots:
# d = dist[2761]
d = dist[79530]
n_nat = np.arange(len(dist[0]))
density_n_nat_d = (3*n_nat[1:])/(4*np.pi*d[1:]**3)
dens_values = np.where(density_n_nat_d <crit_density)[0]


fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(1,1,1)
# ax.scatter(np.arange(len(density_n_nat_d)),density_n_nat_d,marker='.')
ax.plot(n_nat[1:],density_n_nat_d)
ax.plot(n_nat[1:],crit_density*np.ones(len(n_nat[1:])))
ax.plot(n_nat[1:],(len(box)/(box.size()**3))*np.ones(len(n_nat[1:])))

# ax.scatter(d[1:],density_n_nat_d)
#plt.savefig("rad_dens.jpg")
plt.show()

## rad
# Ideal Rad
ideal_rad = (3*n_nat/(4*np.pi*crit_density))**(1/3)


fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(1,1,1)
ax.plot(n_nat[:20],d[:20])
ax.plot(n_nat[:20],ideal_rad[:20])

plt.savefig("rad_dens.jpg")
plt.show()












density_n_nat_d = (3*len(tracers_new[79530]))/(4*np.pi*len(tracers_new[79530])**3)










fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(1,1,1)
# ax.scatter(np.arange(len(density_n_nat_d)),density_n_nat_d,marker='.')
ax.plot(n_nat[1:],density_values[79530])
ax.plot(n_nat[1:],crit_density*np.ones(len(n_nat[1:])))
ax.plot(n_nat[1:],(len(box)/(box.size()**3))*np.ones(len(n_nat[1:])))
ax.axvline(x=len(tracers_new[79530]))
# ax.scatter(d[1:],density_n_nat_d)
#plt.savefig("rad_dens.jpg")
plt.show()