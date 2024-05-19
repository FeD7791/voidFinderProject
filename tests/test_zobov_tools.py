import numpy as np
import scipy as sp
from voidfindertk import box as b
from voidfindertk.zobovvf.zobov_tools import ZobovVoids,ZobovVF ,calculate_tracers_inside_void
from voidfindertk import data_box
##Run just one pytest tests/test_zobov_tools.py -k 'test_calculate_tracers_inside_void'
##Run just one pytest tests/test_zobov_tools.py -k 'test_ZobovVF'

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

def test_ZobovVF(mkbox, make_zobov_voids_params):
    model = ZobovVF()
    box = mkbox(seed=42, size=1000) #Make a box
    dbox = data_box.DataBox(box)
    params = make_zobov_voids_params(n_voids=100,CoreParticle=len(dbox.box)) #Create params for zobov
    
    #model_find step
    voids = ZobovVoids(**params) #Create voids
    voids._tracers_in_void = calculate_tracers_inside_void(dbox.box,voids)
    voids = {'voids':voids}
    ### End model_find step
    #Test mk_vbox
    assert model.mk_vbox(voids,dbox).shape == (len(dbox.box), len(voids['voids'])) 
    #Test get_void_mass
    #sum(voids._tracers_in_void,[]) transform list of list into one list
    index = sum(voids['voids']._tracers_in_void,[])
    a = sum(model.get_void_mass(voids=voids['voids'],llbox=dbox))
    b = sum(box.m.value[index])
    #Round result cause some precision is lost with dataframe
    assert round(a,6) == round(b,6)