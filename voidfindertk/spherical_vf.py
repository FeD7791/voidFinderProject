from abc import ABC
import os
import ctypes
import numpy as np
from . import spherical_voids
from . import spherical_interface
from . import spherical_classifier
class SphericalVF(ABC):
    def __init__(self):
        pass

    def preprocess(databox):
        return databox

    def model_find(llbox):
        sp_void = spherical_interface(llbox.box)
        return {'voids':sp_void}
    
    def mk_vbox(databox,voids,llbox):
        voids = voids['voids']
        databox.__dict__.pop('_len')
        voids.__dict__.pop('_void_len')
        sparse_m = spherical_classifier(**databox.__dict__, **voids.__dict__)
        return sparse_m
