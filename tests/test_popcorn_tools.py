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