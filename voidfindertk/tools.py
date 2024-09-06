import grispy as gsp
import numpy as np
import matplotlib.pyplot as plt
import pathlib
import pandas as pd
from voidfindertk import io

def compute_void_size_function(delta, rsph1, box):
    """
    Computes the void size function

    Parameters
    ----------
        delta : float
        Integrated density Contrast 
        rsph1 : array
        Array of radius of each void in Mpc
        n_tracers : int
        Number of tracers in the whole box.
        box_size : int
        Length in Mpc of the box
    """
    pi = np.pi
    # Volume of the simulation
    # vol = 1000**3
    vol = box.size()**3

    # Mean density of tracers
    rhomed = 6784818/ vol

    # N sequence
    N = np.concatenate([np.arange(1, 20,2),np.arange(20,30,4)])# np.arange(6, 11, 2), np.arange(12, 53, 10)])

    # Scaling calculation
    scl = np.log10((3 / (4 * pi) * N / rhomed / (1 + delta))**(1/3))

    mxlg = np.log10(max(rsph1))
    br = np.concatenate([scl[:-1], np.linspace(max(scl), mxlg, 25)])
    
    # Histogram calculation
    h, bin_edges = np.histogram(np.log10(rsph1), bins=br)

    # Compute dlogr
    dlogr = np.diff(br)
    mids = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    density = h / dlogr / vol

    # Remove zeros
    index = np.where(density>0.0)[0] #Non zero elements index

    return mids[index], density[index]

def calculate_r_eff(*,centers,n_neighbors=100,box,delta=-0.9):
    """
    Calculates the radius for a void when the center and box of tracers are
    provided.

    The final radius is the arithmetic mean between density(n,r_n) and
    density(n+1,r_n+1), where density(n,r_n) < crit_density < density(n+1,r_n+1)
    where density(n,r_n) = n/(4*pi*r_n**3)/3 (Spherical Volume)

    Parameters
    ----------
        centers : array of (x,y,z)
        Array of tracers coordinates.

        n_neighbors : int
        Maximun number of tracers in each void used to perform the search.

        box : Box Object
        Object with the tracer properties.

        delta : float -1< delta < 0
        Integrated density contrast.

    Returns
    -------
        rad_new : numpy array
        Radius calculated for each void.

        tracers_new : tuple of lists
        List of indexes listing the tracers inside each void

        density_values : list of arrays
        For each void gives 1 < n <= n_neighbors density(n,r_n) calculations.
    """
    xyz = np.column_stack((box.x.value,box.y.value,box.z.value))
    grid = gsp.GriSPy(xyz)
    # For each center, get the distance for the n nearest tracers and their index

    dist,nn=grid.nearest_neighbors(centres=centers,n=n_neighbors)

    tracers_new = []
    rad_new = []
    density_values = []
    n_nat = np.arange(1,n_neighbors+1)
    crit_density = (1+delta)*(len(box)/(box.size()**3))
    for n,d in enumerate(dist):
        # Find density values for n_nat particles at radius d
        density_n_nat_d = (3*n_nat)/(4*np.pi*d**3)
        density_values.append(density_n_nat_d)
        # Find all density values that are less than crit_density
        dens_values = np.where(density_n_nat_d <crit_density)[0]
        if(len(dens_values)==0):# This means that all calculated densities are above crit density, probably not a void
            pass
        else:
            # From the values that fulfill the latter condition find the index of the
            # value with max radii
            dist_max_index = np.where(d == max(d[dens_values]))[0][0]


            # Final radii is half distance between distk_max_index and dist_max_index +1
            try:
                rad_new.append((d[dist_max_index+1]+d[dist_max_index])/2)
            except IndexError:
                print(
                    f"All density values under crit density for void center {n} increase n_neighbors to find right number of tracers in void"
                    )
            else:
            # Give the all the tracers whose distance is lesser than the found radii
                tracers_new.append(nn[n][:dist_max_index])
    ##For testing
    #Demostrar que la densidad , calculada con rad_new y len(tracers_new) da menor a la densidad critica
    #Demostrar que cada trazador generado esta dentro del radio dado, lo que significa que se respeta el orden del box
    #Para finders como el esferico, demostrar que se obtienen los mismos trazadores
    #Para finders como el esferico, demostrar que la cantidad de voids son los mismos
    return rad_new,tuple(tracers_new),density_values

