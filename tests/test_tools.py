import ctypes
import numpy as np
import pytest
from voidfindertk.tools import ctypes_input_params_builder, ctypes_output_builder, process_output_from_finder, preprocess_data_box
from voidfindertk import box
from voidfindertk.models import DataBox

def test_ctypes_input_params_builder():
    params_list_double = [
        "RadIncrement", "BoxSize",
        "MaxRadiusSearch", "ProxyGridSize",
        "FracRadius", "DeltaThreshold",
        "DeltaSeed", "OverlapTol", 
        "Redshift", "OmegaMatter",
        "OmegaLambda", "Hubble",
        "FidOmegaMatter", "FidOmegaLambda",
        "FidHubble", "MinProfileDist",
        "MaxProfileDist", "ScalePos",
        "ScaleVel", "InnerShellVel",
        "OuterShellVel"
    ]
    params_list_int = [
        "FormatTracers", "NumFiles",
        "NumRanWalk", "OMPcores",
        "RSDist", "GDist",
        "WriteProfiles", "NumProfileBins"
    ]

    double_values = len(params_list_double)*[ctypes.c_double]
    int_values = len(params_list_int)*[ctypes.c_int]
    params = params_list_double + params_list_int
    values = len(params)*[0]
    ctypes_values = double_values + int_values
    params_dict = dict(list(zip(params,values)))
    output = ctypes_input_params_builder('spherical',**params_dict)
    assert output['params_dict'] == params_dict
    assert output['InputParams_class']._fields_ == list(zip(params,ctypes_values)) 


def test_ctypes_output_builder(make_spherical_voids_params):
    params = ["n_voids", "Rad", "Pos", "Vel", "Dtype", "Delta", "Poisson", "Nran"]
    ctypes_values = [ctypes.c_int] + [ctypes.c_float] + 2*[ctypes.c_float*3] + 3*[ctypes.c_float] + [ctypes.c_int]
    voids_output = list(zip(params,ctypes_values))
    output = ctypes_output_builder('spherical')
    assert output['voidArray_class']._fields_[0][1]._type_._fields_ == voids_output   

def test_process_output_from_finder():
    array_of_voids = [
        [
            [
                [
                    2,
                    'rad',
                    ['x','y','z'],
                    ['vx','vy','vz'],
                    'dtype',
                    'delta',
                    'poisson',
                    'nran'
                    ]]]]
    output = process_output_from_finder('spherical',np.array(array_of_voids, dtype='object'))
    assert output['rad'][0] == 'rad'
    assert output['x_void'][0] == 'x'
    assert output['y_void'][0] == 'y'
    assert output['z_void'][0] == 'z'
    assert output['vel_x_void'][0] == 'vx'
    assert output['vel_y_void'][0] == 'vy'
    assert output['vel_z_void'][0] == 'vz'
    assert output['delta'][0] == 'delta'
    assert output['dtype'][0] == 'dtype'
    assert output['poisson'][0] == 'poisson'
    assert output['nran'][0] == 'nran'

def test_preprocess_data_box():
    input_data = {
        'x':np.arange(0,1000,1),
        'y':np.arange(0,1000,1),
        'z':np.arange(0,1000,1),
        'vx':np.arange(0,1000,1),
        'vy':np.arange(0,1000,1),
        'vz':np.arange(0,1000,1),
        'm':np.arange(0,1000,1),
                  }
    b = box.Box(**input_data) # creating box from input data
    d_b = data_box.DataBox(b) # creating data_box from box
    # Test function without keywords
    output = preprocess_data_box(d_b) # m should be bigger than 50, m > 50
    assert len(output.box.m) == len(input_data['m'])
    assert len(output.box.x) == len(input_data['x'])
    assert len(output.box.y) == len(input_data['y'])
    assert len(output.box.z) == len(input_data['z'])
    assert len(output.box.vx) == len(input_data['vx'])
    assert len(output.box.vy) == len(input_data['vy'])
    assert len(output.box.vz) == len(input_data['vz'])
    # Test min mass
    output = preprocess_data_box(d_b,m_min = 50) # m should be bigger than 50, m > 50
    assert len(output.box.m) == len(input_data['m'][input_data['m'] >= 50])
    assert len(output.box.x) == len(input_data['x'][input_data['x'] >= 50])
    assert len(output.box.y) == len(input_data['y'][input_data['y'] >= 50])
    assert len(output.box.z) == len(input_data['z'][input_data['z'] >= 50])
    assert len(output.box.vx) == len(input_data['vx'][input_data['vx'] >= 50])
    assert len(output.box.vy) == len(input_data['vy'][input_data['vy'] >= 50])
    assert len(output.box.vz) == len(input_data['vz'][input_data['vz'] >= 50])

    # Test max mass
    output = preprocess_data_box(d_b,m_max = 50) 
    assert len(output.box.m) == len(input_data['m'][input_data['m'] <= 50])
    assert len(output.box.x) == len(input_data['x'][input_data['x'] <= 50])
    assert len(output.box.y) == len(input_data['y'][input_data['y'] <= 50])
    assert len(output.box.z) == len(input_data['z'][input_data['z'] <= 50])
    assert len(output.box.vx) == len(input_data['vx'][input_data['vx'] <= 50])
    assert len(output.box.vy) == len(input_data['vy'][input_data['vy'] <= 50])
    assert len(output.box.vz) == len(input_data['vz'][input_data['vz'] <= 50])

    # Test m_max bigger than len(output.box.m)
    output = preprocess_data_box(d_b,m_max = 2500)
    assert len(output.box.m) == len(input_data['m'])

    # test not duplicates
    input_data = {
        'x':np.array([2,2,2,3,3,3]),
        'y':np.array([2,2,2,3,3,3]),
        'z':np.array([2,2,2,3,3,3]),
        'vx':np.array([2,2,2,3,3,3]),
        'vy':np.array([2,2,2,3,3,3]),
        'vz':np.array([2,2,2,3,3,3]),
        'm':np.array([2,2,2,3,3,3]),
                  }
    b = box.Box(**input_data) # creating box from input data
    d_b = data_box.DataBox(b) # creating data_box from box
    output = preprocess_data_box(d_b) 
    assert all(output.box.m.value == np.array([2,3]))
    assert all(output.box.x.value == np.array([2,3]))
    assert all(output.box.y.value == np.array([2,3]))
    assert all(output.box.z.value == np.array([2,3]))
    assert all(output.box.vx.value == np.array([2,3]))
    assert all(output.box.vy.value == np.array([2,3]))
    assert all(output.box.vz.value == np.array([2,3]))
    