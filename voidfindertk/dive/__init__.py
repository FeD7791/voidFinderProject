"""DIVE: Detection, Identificacion and Voids Examination

DIVE is a reimplementation of VIDE void detection algoritm in pure Python on
top of the classic ZOBOV code.

"""

import datetime as dt
import os
import pathlib
import shutil
import tempfile

import numpy as np

from . import _postprocessing
from . import _wrapper as _wrap
from ..zobov import ZobovVF


class DiveVF(ZobovVF):

    def __init__(self, *, dive_coso, **kwargs):
        super().__init__(**kwargs)
        self._dive_coso = dive_coso

    @property
    def dive_coso(self):
        return self._dive_coso

    def build_voids(self, model_find_parameters):
        particle_by_voids, extra = super().build_voids(model_find_parameters)
