#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Rava Jorge Federico, Gualpa Sebastian
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
import numpy as np

from astropy import units as u

import pytest

from voidfindertk.box import Box


def test_Box_initialization(mkbox):
    box = mkbox(seed=42, size=1000)

    assert box.x.unit == u.Mpc
    assert box.y.unit == u.Mpc
    assert box.z.unit == u.Mpc
    assert box.vx.unit == u.Mpc / u.h
    assert box.vy.unit == u.Mpc / u.h
    assert box.vz.unit == u.Mpc / u.h
    assert box.m.unit == u.M_sun
    assert len(box) == 1000
    assert repr(box) == "<Box size=1000>"
    assert isinstance(box.x, np.ndarray)
    assert isinstance(box.y, np.ndarray)
    assert isinstance(box.z, np.ndarray)
    assert isinstance(box.vx, np.ndarray)
    assert isinstance(box.vy, np.ndarray)
    assert isinstance(box.vz, np.ndarray)
    assert isinstance(box.m, np.ndarray)


def test_Box_different_length_for_some_dimension(mkbox_params):
    params = mkbox_params(seed=42, size=1000)
    params["x"] = params["x"][:-1]

    with pytest.raises(ValueError, match="Arrays should be of the same size"):
        Box(**params)


def test_box_equality(mkbox):
    box1 = mkbox(seed=42, size=1000)
    box2 = mkbox(seed=42, size=1000)
    box3 = mkbox(seed=43, size=1000)
    box4 = mkbox(seed=42, size=999)
    assert box1 == box2
    assert not box1 == box3
    assert not box4 == box1
