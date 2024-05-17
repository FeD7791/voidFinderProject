import numpy as np
import scipy as sp
from voidfindertk.sphericalvf import spherical_tools
from voidfindertk import comparison_tools


## pytest tests/test_comparison_tools.py -k 'test_find_radius'
def test_find_radius(mkbox,make_spherical_voids_params):
    #box
    box = mkbox(seed=42, size=1000)
    box_0 = np.array([box.x[0].value,box.y[0].value,box.z[0].value])
    box_50 = np.array([box.x[50].value,box.y[50].value,box.z[50].value]) 
    #sparse matrix 
    #this part represent the whole box against 2 voids found
    sparse_matrix = np.zeros((len(box), 2))# the matrix only have two ones in first column position 0 and 50
    sparse_matrix[0,0] = 1
    sparse_matrix[50,0] = 1
    #voids
    params = make_spherical_voids_params(n_voids=2)
    spherical_voids = spherical_tools.SphericalVoids(**params)
    void_center = np.array([
        spherical_voids.x_void[0].value,
        spherical_voids.y_void[0].value,
        spherical_voids.z_void[0].value,
        ])
    max_distance = sp.spatial.distance.pdist([box_0,box_50,void_center]).sort(-1)
    assert comparison_tools.find_radius(
        box=box,
        voids=spherical_voids,
        sparse=sparse_matrix)[0] == max_distance[-1]

    