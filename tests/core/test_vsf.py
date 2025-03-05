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

"""Tests for the voidfindertk.core.vsf."""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np

import pytest

from voidfindertk.core import vsf

# =============================================================================
# TESTS
# =============================================================================

rng = np.random.default_rng(seed=50)


@pytest.mark.parametrize(
    "radius, delta",
    [
        # Uniform Distribution
        (rng.uniform(10, 100, size=(1000, 1)), -0.7),
        (rng.uniform(10, 100, size=(1000, 1)), -0.8),
        (rng.uniform(10, 100, size=(1000, 1)), -0.9),
        # Normal Distribution
        (rng.normal(loc=30, scale=15, size=(1000, 1)), -0.7),
        (rng.normal(loc=30, scale=15, size=(1000, 1)), -0.8),
        (rng.normal(loc=30, scale=15, size=(1000, 1)), -0.9),
    ],
)
def test_void_size_function(mkbox, radius, delta):
    box = mkbox(seed=42)
    log_of_radius, counts, delta_out = vsf.void_size_function(
        radius=np.ravel(radius), delta=delta, box=box
    )
    assert delta_out == delta
    assert len(log_of_radius) > 0, "log_of_radius should not be empty"
    assert len(counts) > 0, "counts should not be empty"
    assert np.all(log_of_radius > 0), "log_of_radius values should be positive"
    assert np.all(counts > 0), "counts should have non-zero values"
