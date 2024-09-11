import numpy as np
from voidfindertk import box
import grispy as gsp

Npoints = 100 ** 3
Ncentres = 10
dim = 3
Lbox = 1000.0

centers = np.array([np.array([10,10,10]),np.array([20,20,20]),np.array([30,30,30])])

rad_max = np.array([20,40])
#build box

def build_cloud(
        *,
        seed=2,
        lmin=0,
        lmax=1000,
        n_points=100**3
        ):
    np.random.seed(seed)
    #Create point cloud
    cloud = np.random.uniform(lmin, lmax, size=(n_points, 3))
    return cloud

kloud = build_spherical_void(delta=-0.7,centers=centers,radii=30,cloud=cloud)

def build_spherical_void(
        delta,
        centers:np.ndarray,
        radii:float,
        cloud,
        ):
    cloud_volume = round(np.max(cloud)-np.min(cloud))**3
    cloud_density = len(cloud)/cloud_volume
    density_voids = (1+delta)*cloud_density

    #Find neares neighbors
    grid = gsp.GriSPy(cloud)
    dist,index = grid.bubble_neighbors(centers, distance_upper_bound=radii)

    #Calculate right number of tracers so the void gets at the desired density
    n = round(density_voids*(4/3)*np.pi*(radii**3))


    new_index = [i[:len(i)-n+1] for i in index] #+1 to ensure we remove more particles so density of the void < dens_voids

    # Remove from cloud the indicated indexes
    mask = np.ones(len(cloud),dtype=bool)
    for i in new_index:
        mask[i] = False

    cloud_with_voids = cloud[mask]
    return cloud_with_voids









#Spherical voids
b_dist,b_index = grid.bubble_neighbors(centers[:2], distance_upper_bound=rad_max)



#Make spherical voids
mask = np.ones(len(data),dtype=bool)
for index in b_index:
    mask[index] = False

xyz = data[mask]

#Spherical with high density


#add density
data = np.random.uniform(300, 305, size=(0, dim))





grid2 = gsp.GriSPy(xyz)
b_dist,b_index = grid2.shell_neighbors(np.array([centers[2]]), distance_upper_bound=30, distance_lower_bound=5)

mask = np.ones(len(xyz),dtype=bool)
for index in b_index:
    mask[index] = False

xyz2 = xyz[mask]

x,y,z = np.hsplit(xyz2,3)
dbox = box.Box(**{
    'x':x,
    'y':y,
    'z':z,
    'vx':x,
    'vy':y,
    'vz':z,
    'm':z
    })

xyz3 = np.vstack((xyz2,data))

x,y,z = np.hsplit(xyz3,3)
dbox = box.Box(**{
    'x':x,
    'y':y,
    'z':z,
    'vx':x,
    'vy':y,
    'vz':z,
    'm':z
    })