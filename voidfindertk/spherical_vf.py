from abc import ABC
import os
import ctypes
import numpy as np
from . import spherical_voids
from . import spherical_interface
from . import spherical_classifier
from .models import ModelABC
class SphericalVF(ModelABC):
    def __init__(self):
        pass

    def preprocess(self, databox):
        return databox

    def model_find(self, llbox):
        sp_void = spherical_interface.spherical_void_finder(llbox.box)
        return {'voids':sp_void}
    
    def mk_vbox(self, databox,voids,llbox):
        voids = voids['voids']
        databox.box.__dict__.pop('_len')
        voids.__dict__.pop('_void_len')
        sparse_m = spherical_classifier.Box_SphericalVoids_Classifier(**databox.box.__dict__, **voids.__dict__)._sparse_matrix(0.0)
        return sparse_m
