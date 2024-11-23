#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023 - 2024, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
#                            Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.

# =============================================================================
# DOCS
# =============================================================================

"""
Test ZOBOV posprocessing modules

"""

# =============================================================================
# IMPORTS
# =============================================================================

import pathlib

import pytest

import numpy as np

import pandas as pd

import sh

from voidfindertk.zobov._zobov import _Files
from voidfindertk.zobov import _zb_postprocessing




# =============================================================================
# TESTS
# =============================================================================



def test_get_tracers_in_zones():
    path = pathlib.Path("/home/jorgefederico/updates/vftk_1109/voidFinderProject/runz/tmpek1e95uj2024-11-21T20:13:55.126949+00:00")
    # Dataframe with out file .txt values
    out = pd.read_csv(path/_Files.OUTPUT_JOZOV_VOIDS_DAT,header=1,delim_whitespace=True)
    # Tracers that are in each zone
    tracers = _zb_postprocessing._get_tracers_in_zones(
        tracers_in_zones_path=path/_Files.PARTICLES_VS_ZONES_ASCII)

    # This have the same order:
    core_particle = np.array([t[0] for t in tracers]) 
    out_core_particle = np.array(out['CoreParticle'])

    indx = np.argsort(out_core_particle)

    # the number of tracers in tracers should be the same as out['Zone#Part']
    # when they are in the right order
    number_of_tracers = np.array(list(map(len,tracers)))

    assert np.all(core_particle==out_core_particle[indx])
    assert np.all(number_of_tracers==np.array(out['Zone#Part'][indx]))

def test_get_zones_in_void():
    path = pathlib.Path("/home/jorgefederico/updates/vftk_1109/voidFinderProject/runz/tmpek1e95uj2024-11-21T20:13:55.126949+00:00")
    # Dataframe with out file .txt values
    out = pd.read_csv(path/_Files.OUTPUT_JOZOV_VOIDS_DAT,header=1,delim_whitespace=True)

    void_number = np.array(out['Void#'])
    zones = _zb_postprocessing._get_zones_in_void(path/_Files.ZONES_VS_VOID_ASCII)
    zone_number = np.array([z[0] for z in zones]) -1
    number_of_zones = np.array(out['Void#Zones'])
    # You can think of 'two' columns 1: #Zone 2: then #Zone + other zones
    number_of_zones_from_module = np.array(list(map(len,zones))) -1
    import ipdb; ipdb.set_trace()
    assert np.all(void_number==zone_number)
    assert np.all(number_of_zones==number_of_zones_from_module)

def test_get_tracers_in_voids(zobov_model_builder):

    extra = zobov_model_builder()
    path = extra["files_directory_path"]


    # path = pathlib.Path("/home/jorgefederico/updates/vftk_1109/voidFinderProject/runz/tmpek1e95uj2024-11-21T20:13:55.126949+00:00")
    out = pd.read_csv(path/_Files.OUTPUT_JOZOV_VOIDS_DAT,header=1,delim_whitespace=True)
    properties_df, tracers_in_void = _zb_postprocessing.get_tracers_in_voids(
        properties_dataframe=out,
        tracers_in_zones_path=path/_Files.PARTICLES_VS_ZONES_ASCII,
        zones_in_void_path=path/_Files.ZONES_VS_VOID_ASCII
        )
    number_tracers_in_void = np.array(list(map(len,tracers_in_void)))
    df_n_tracers_in_void = np.array(properties_df['Void#Part'])
    sh.rm('-rf',path)
    assert np.all(number_tracers_in_void==df_n_tracers_in_void)