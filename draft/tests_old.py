# test_comparison_tools.py
import numpy as np
import scipy as sp
from voidfindertk.sphericalvf import spherical_tools
from voidfindertk import tools


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
    assert tools.find_radius(
        box=box,
        voids=spherical_voids,
        sparse=sparse_matrix)[0] == max_distance[-1]

# test pipeline
from voidfindertk.models import DataBox
from voidfindertk.sphericalvf.spherical_tools import SphericalVF, spherical_void_finder, SphericalVoids


def test_sphericalvf_class(mkbox, make_spherical_voids_params):
    box = mkbox(seed=42, size=1000) #Build artificial box of size 1000
    params = make_spherical_voids_params() #Build params for spherical void output
    data_box = DataBox(box) #Transform box to universal Box
    spherical_void = SphericalVoids(**params) #Build artificial void output using params
    model_find_output = {'voids':spherical_void} #Artificial output for mk_vbox step
    model = SphericalVF()
    output = model.mk_vbox(model_find_output,data_box) #Build artificial mk_vbox
    assert box == data_box.box #For spherical, data_box is the same as box 
    assert data_box == model.preprocess(databox=data_box) #preprocess in spherical won't modify data_box 
    assert output['box'] == data_box.box 
    assert output['voids'] == spherical_void
    assert (len(box),len(spherical_void)) == output['sparse'].shape

#popcorn tools
import os
from voidfindertk.popcornvf.popcorn_tools import popcorn_svf_input_data_builder, config_file_maker
from voidfindertk.io import read_table
import numpy as np

##Run just one pytest tests/test_popcorn_tools.py -k 'test_config_file_maker'

def test_popcorn_svf_input_data_builder(mkbox):
    box = mkbox(seed=3, size=1000) #Build artificial box of size 1000
    path = os.path.dirname(os.path.realpath(__file__)) #Path of the thest file
    path_vftk = os.path.split(path) #path_vftk[0] path of the voidfinderproject folder
    path_file = os.path.join(
        path_vftk[0],
        'voidfindertk',
        'popcorn',
        'popcorn_void_finder',
        'Source') #Go from test to source
    popcorn_svf_input_data_builder(box=box)
    path_input_data= os.path.join(path_file,'box.dat') #Path where the file is
    box2 = read_table(
        path_input_data, 
        float_precision='round_trip' # Pandas uses a dedicated decimal-to-binary converter that sacrifices perfect accuracy for the sake of speed
        ) #Read created data, and put it in a box2
    os.remove(path=path_input_data) #Remove the created file
    assert np.array_equal(box.x[0].value,box2.x[0].value)
    assert all(box.x.value == box2.x.value)
    assert all(box.y.value == box2.y.value)
    assert all(box.z.value == box2.z.value) 
    assert all(box.vx.value == box2.vx.value) 
    assert all(box.vy.value == box2.vy.value)
    assert all(box.vz.value == box2.vz.value)
    assert all(box.m.value == box2.m.value) 

def test_config_file_maker():
    params = {
        'TRSFILE':'/home/some/path',
        'FILEFMT':'ASCII',
        'NUM_FILE': 1,
        'POPFILE' : '../MWE/clean_popvds.dat',
        'SPHFILE' : '../MWE/sphvds.dat',
        'AUXFILES' : 'true',
        'RAWPOPFILE' : '../MWE/raw_popvds.dat',
        'RAWSPHFILE' : '../MWE/raw_sphv.dat',
        'PAIRSFILE' : '../MWE/pairs_raw_popvds.dat',
        'BOXSIZE' : 1000.0,
        'DENSTH' : -0.9,
        'MINRADIUS' : 5,
        'MAXRADIUS' : 100,
        'EPS' : 0.00001,
        'MASSMIN' : 50, 
    }
    path = os.path.dirname(os.path.realpath(__file__)) 
    params_file = [f'{key} = {value}' for key,value in params.items()] #compare with output of file
    config_file_maker(params, path) #this will create the file
    with open(os.path.join(path,'vars.conf'),'r') as vars:
        output = vars.readlines(-1)
    output.pop(0) #removing [INPUT_PARAMS]
    output_clean = [output[i].replace('\n','') for i in range(15)]
    os.remove(path=os.path.join(path,'vars.conf')) #Remove the created file
    assert output_clean == params_file


#shperical tools
from voidfindertk.sphericalvf import spherical_tools
from astropy import units as u
import pytest
from voidfindertk.sphericalvf.spherical_tools import SphericalVoids
import numpy as np
import pandas as pd

def test_SphericalVoids(make_spherical_voids_params):
    params = make_spherical_voids_params(n_voids=500)
    spherical_voids = SphericalVoids(**params)
    assert spherical_voids.rad.unit == u.Mpc
    assert spherical_voids.x_void.unit == u.Mpc
    assert spherical_voids.y_void.unit == u.Mpc
    assert spherical_voids.z_void.unit == u.Mpc
    assert spherical_voids.vel_x_void.unit == u.Mpc/u.h
    assert spherical_voids.vel_y_void.unit == u.Mpc/u.h
    assert spherical_voids.vel_z_void.unit == u.Mpc/u.h
    assert len(spherical_voids) == 500
    assert str(spherical_voids) == "<SphericalVoids size=500>"
    #wont test slice method because i dont know if i will keep it

# test wrapper

from voidfindertk.zobov._wrapper import run_vozinit


def test_run_vozinit():

    input_params = {
        'position_file': 'tracers_zobov.raw',
        'buffer_size': '0.8',
        'box_size': '500',
        'number_of_divisions': '2',
        'suffix_describing_this_run': 'output'
    }
    len_box = '700070'
    generic_output_run_vozinit = [
        f'boz = {input_params["box_size"]}',
        f'np = {len_box}',
        f'np: {len_box}, x: 1.06e-06,0.999999; y: 8e-07,1; z: 8.4e-07,1',
        's = 0.0119048, bf = 0.8, g = 0.799911.',
        'b=(0,0,0), c=(0.25,0.25,0.25), nvp=85789, nvpbuf=700070',
        'b=(0,0,1), c=(0.25,0.25,0.75), nvp=82254, nvpbuf=700070',
        'b=(0,1,0), c=(0.25,0.75,0.25), nvp=92853, nvpbuf=700070',
        'b=(0,1,1), c=(0.25,0.75,0.75), nvp=83614, nvpbuf=700070',
        'b=(1,0,0), c=(0.75,0.25,0.25), nvp=88198, nvpbuf=700070',
        'b=(1,0,1), c=(0.75,0.25,0.75), nvp=84321, nvpbuf=700070',
        'b=(1,1,0), c=(0.75,0.75,0.25), nvp=96789, nvpbuf=700070',
        'b=(1,1,1), c=(0.75,0.75,0.75), nvp=86253, nvpbuf=700070',
        'Nvp range: 82254,96789',
        'Nvpbuf range: 700070,700070',
        'Writing script file to scroutput.',
        ''
    ]

# zobov tools
import numpy as np
import scipy as sp
from voidfindertk import box as b
from voidfindertk.zobovvf.zobov_tools import ZobovVoids,ZobovVF ,calculate_tracers_inside_void
from voidfindertk.models import DataBox
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
    dbox = DataBox(box)
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