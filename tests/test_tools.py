import ctypes
import numpy as np
import pytest
from voidfindertk.tools import ctypes_input_params_builder, ctypes_output_builder

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