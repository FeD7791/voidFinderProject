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

"""Test for voidfindertk.core.box."""

# =============================================================================
# IMPORTS
# =============================================================================

from astropy import units as u

import numpy as np

import pytest

from voidfindertk.core.box import Box

# =============================================================================
# TESTS
# =============================================================================


def test_box_initialization(mkbox):
    box = mkbox(seed=42, size=1000)
    assert box.x.unit == u.Mpc
    assert box.y.unit == u.Mpc
    assert box.z.unit == u.Mpc
    assert box.vx.unit == u.Mpc / u.second
    assert box.vy.unit == u.Mpc / u.second
    assert box.vz.unit == u.Mpc / u.second

    assert len(box) == 1000
    assert repr(box) == "<Box size=1000>"
    assert isinstance(box.x, np.ndarray)
    assert isinstance(box.y, np.ndarray)
    assert isinstance(box.z, np.ndarray)
    assert isinstance(box.vx, np.ndarray)
    assert isinstance(box.vy, np.ndarray)
    assert isinstance(box.vz, np.ndarray)
    assert isinstance(box.m, np.ndarray)


def test_box_different_length_tracers(mkbox_params):
    params = mkbox_params(seed=42, size=1000)
    params["x"] = params["x"][:-1]

    with pytest.raises(ValueError, match="Arrays should be of the same size"):
        Box(**params)


def test_box_equality(mkbox):
    box1 = mkbox(seed=42, size=1000)
    box2 = mkbox(seed=42, size=1000)
    box3 = mkbox(seed=42, size=999)
    assert box1 == box2
    assert not box1 == box3


def test_box_mass_cutoff(mkbox):
    m_threshold = 10
    box = mkbox(seed=42, size=1000, mass_scale=15)
    new_box = box.mass_cutoff(mass_threshold=m_threshold)

    assert m_threshold <= np.min(new_box.m)
    assert len(new_box) < len(box)
    assert new_box.max_ <= box.max_
    assert box.min_ <= box.max_
