#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Rava Jorge Federico, Gualpa Sebastian
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
import io

import numpy as np

import pandas as pd

import pytest

from voidfindertk.box import Box



@pytest.fixture(scope="session")
def mkbox_params(): #aca vas a retornar la funcion maker que la podes modificar
    def _maker(
        *,
        seed=None,
        coordinates_scale=500,
        velocity_scale=220,
        mass_scale=12,
        size=1000,
    ):
        rng = np.random.default_rng(seed=seed)
        params = {
            "x": rng.uniform(0,coordinates_scale,size=size),
            "y": rng.uniform(0,coordinates_scale,size=size),
            "z": rng.uniform(0,coordinates_scale,size=size),
            "vx": rng.uniform(0,velocity_scale,size=size),
            "vy": rng.uniform(0,velocity_scale,size=size),
            "vz": rng.uniform(0,velocity_scale,size=size),
            "m": rng.uniform(0,mass_scale,size=size)
        }
        return params

    return _maker


@pytest.fixture(scope="session")
def mkbox(mkbox_params):
    def _maker(**kwargs):
        params = mkbox_params(**kwargs)
        return Box(**params)

    return _maker


@pytest.fixture(scope="session")
def random_buffer():
    def _maker(*, empty_row=False, sp_characters=None):
        data = np.random.random((100, 7))
        df = pd.DataFrame(data)
        n = 0
        if empty_row:
            df.loc[100] = [np.nan] * 7
            n = 1
        elif sp_characters:
            for i in np.arange(len(sp_characters)):
                df.loc[100 + n + i] = [sp_characters[i]] * 7
        src = df.to_csv(
            sep=" ",
            index=False,
            float_format="%.5f",
            header=False,
            na_rep=np.nan,
        )
        buff = io.StringIO(src)
        return buff

    return _maker

@pytest.fixture
def make_spherical_voids_params():
    #aca vas a retornar un valor, no una funcion
    def _maker(**kwargs):
        params = {
            'n_voids':1000,
            'rad_scale': 15,
            'xyz_void_max_scale':500,
            'vel_xyz_void_max_scale':200,
            'min_delta':-0.95,
            'max_delta':-0.90,
            'min_poisson':-0.5,
            'max_poisson':0.5,
            'dtype': 0.5,
            'nran': 200,
            'seed':42
        }
        for key,value in kwargs.items():
            params[key] = value

        rng = np.random.default_rng(seed=params['seed'])
        void_params = {
            'rad' : rng.uniform(0,params['rad_scale'],params['n_voids']),
            'x_void' : rng.uniform(0,params['xyz_void_max_scale'],params['n_voids']),
            'y_void' : rng.uniform(0,params['xyz_void_max_scale'],params['n_voids']),
            'z_void' : rng.uniform(0,params['xyz_void_max_scale'],params['n_voids']),
            'vel_x_void' : rng.uniform(0,params['vel_xyz_void_max_scale'],params['n_voids']),
            'vel_y_void' : rng.uniform(0,params['vel_xyz_void_max_scale'],params['n_voids']),
            'vel_z_void' : rng.uniform(0,params['vel_xyz_void_max_scale'],params['n_voids']),
            'delta' : rng.uniform(params['min_delta'],params['max_delta'],params['n_voids']),
            'poisson' : rng.uniform(params['min_poisson'], params['max_poisson'],params['n_voids']),
            'dtype' :  rng.uniform(0, params['dtype'],params['n_voids']),
            'nran' : rng.uniform(0, params['nran'],params['n_voids'])
        }
        return void_params
    
    return _maker

