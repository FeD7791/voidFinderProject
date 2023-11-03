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
def mkbox_params():
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
            "x": coordinates_scale * rng.random(size=size),
            "y": coordinates_scale * rng.random(size=size),
            "z": coordinates_scale * rng.random(size=size),
            "vx": velocity_scale * rng.random(size=size),
            "vy": velocity_scale * rng.random(size=size),
            "vz": velocity_scale * rng.random(size=size),
            "m": mass_scale * rng.random(size=size),
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
            df.loc[100] = [np.nan]*7
            n = 1
        elif sp_characters:
            for i in np.arange(len(sp_characters)):
                df.loc[100 + n + i] = [sp_characters[i]]*7
        src = df.to_csv(sep=" ", index=False, float_format="%.5f", header=False, na_rep=np.nan)
        buff = io.StringIO(src)
        return buff
    
    return _maker
