import numpy as np
import scipy as sp
from voidfindertk import box as b
from voidfindertk.zobovvf.zobov_tools import ZobovVoids, calculate_tracers_inside_void

##Run just one pytest tests/test_zobov_tools.py -k 'test_calculate_tracers_inside_void'

def test_calculate_tracers_inside_void(mkbox, make_zobov_voids_params):
    box = mkbox(seed=42, size=1000) #Make a box
    params = make_zobov_voids_params(n_voids=2) #Create params for zobov
    voids = ZobovVoids(**params) #Create voids
    for i in range(1000): #Fill with artificial positions of the tracers
        box.x.value[i] = i 
        box.y.value[i] = i 
        box.z.value[i] = i 
    #Change values for the void so the centers are in the 0 and 999 tracers
    voids.CoreParticle[0] = 0
    voids.CoreParticle[1] = 999
    # The number of particles in voids 0 and 999 are 5 and 3 respectively
    voids.Void_number_Part[0] = 5
    voids.Void_number_Part[1] = 3
    #5 Closest particle of tracer 0 are [1,2,3,4,5] y closest 3 to tracer 9 are [8,7,6] 
    assert calculate_tracers_inside_void(box,voids) == [[1,2,3,4,5],[998,997,996]]
