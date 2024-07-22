import numpy as np

#--------------------------------------------------------------------------
# Spherical reader

def read_sph(filename):
    """
    Reads a data table of spherical voids from an ASCII file (without header) and structures it as a dictionary.

    Parameters:
    - filename (str): the path to the data file.

    Returns:
    - sphvoids (dict): a dictionary with keys corresponding to spherical voids attributes.

    Spherical voids attribute:
    * id: void ID (integer)
    * r: radii of the member spheres (float)
    * x: x-coordinate of the centres of the member spheres (float)
    * y: y-coordinate (float)
    * z: z-coordinate (float)
    * delta: integrated density contrast at void radius (float)
    """

    file = np.loadtxt(filename)

    nsph = len(file)
    print(f"Number of spheres: {nsph}")

    columns = ["id", "r", "x", "y", "z", "delta"]
    sphvoids = {col: file[:,i] for i, col in enumerate(columns)}
    
    return sphvoids

#--------------------------------------------------------------------------
# Popcorn reader

def read_pop(filename):
    """
    Reads a data table of popcorn voids from an ASCII file (without header) and structures it as a dictionary.

    Parameters:
    - filename (str): The path to the data file.

    Returns:
    - popcorn (dict): a dictionary with keys corresponding to popcorn voids attributes.

    Popcorn voids attribute:
    * id: popcorn ID (integer)
    * nmem: number of member spheres (integer)
    * vol: popcorn volume (float)
    * reff: effective radius (float)
    * npart: number of particles inside (integer)
    * flag: for internal control, you can ignore it... (integer)
    * pop: popcorn object with the attributes of the member spheres (dictionary, length: nmem)
        * x: x-coordinate of the centres of the member spheres (float)
        * y: y-coordinate (float)
        * z: z-coordinate (float)
        * r: radii of the member spheres (float)
        * fvol: volume contribution to the total volume (float)
        * level: hierarchy level: 0 for main sphere, 1 for secondary spheres, ... (integer)
    """

    with open(filename, 'r') as file:
        npop = int(file.readline().strip()) # number of popcorn voids (integer)
        print(f"Number of popcorns: {npop}")
        
        popcorn = {'id':[], 'nmem':[], 'vol':[], 'reff':[], 'npart':[], 'flag':[], 'pop':[]}

        for p in range(npop):
            head_pop = list(file.readline().strip().split()) # this line is the header
            id    = int(head_pop[0])
            nmem  = int(head_pop[1])
            vol   = float(head_pop[2])
            npart = int(head_pop[3])
            flag  = int(head_pop[4])

            reff = (vol*3/(4*np.pi))**(1./3.)
	    # NOTE: this is a derived result, not contained in the original catalogue

            # Building the `pop` object:
            x=[]; y=[]; z=[]; r=[]; fvol=[]; level=[]
            if nmem > 0:
                for _ in range(nmem):
                    popmem = list(file.readline().strip().split()) # attributes of the member spheres
                    x.append(float(popmem[0]))
                    y.append(float(popmem[1]))
                    z.append(float(popmem[2]))
                    r.append(float(popmem[3]))
                    fvol.append(float(popmem[4]))
                    level.append(int(popmem[5]))
            pop = {'x':x, 'y':y, 'z':z, 'r':r, 'fvol':fvol, 'level':level} # this is an individual popcorn

            if npart > 0:
                for _ in range(npart):
                    file.readline() # reading inner particles ID (ignoring)

            popcorn['id'].append(id)
            popcorn['nmem'].append(nmem)
            popcorn['vol'].append(vol)
            popcorn['reff'].append(reff)
            popcorn['npart'].append(npart)
            popcorn['flag'].append(flag)
            popcorn['pop'].append(pop)

    return popcorn
