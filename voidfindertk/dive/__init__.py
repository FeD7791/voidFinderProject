#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the
#   Void-Finder-Toolkit Project (https://github.com/FeD7791/voidFinderProject).
# Copyright (c) 2023-2024, Bustillos, J.F.; Cabral, Juan
# License: BSD 3-Clause
# Full Text:
#   https://raw.githubusercontent.com/FeD7791/voidFinderProject/dev/LICENSE.txt

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
