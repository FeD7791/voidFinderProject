#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Rava Jorge Federico, Gualpa Sebastian
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
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


def test_Box_different_length_for_some_dimension(mkbox_params):
    params = mkbox_params(seed=42, size=1000)
    params["x"] = params["x"][:-1]

    with pytest.raises(ValueError, match="Arrays should be of the same size"):
        Box(**params)
