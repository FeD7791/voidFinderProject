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

""" tests for the voidfindertk.core.vfinder_abc

"""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np

import pytest

from voidfindertk.core.vfinder_abc import VoidFinderABC

# =============================================================================
# TESTS
# =============================================================================


def test_VFinderABC(mkbox):
    box = mkbox(seed=42)

    class TestVFinder(VoidFinderABC):
        def preprocess(self, box):
            return [1]

        def model_find(self, preprocess_parameters):
            return preprocess_parameters + [2]

        def build_voids(self, model_find_parameters):
            tracers = (1, 2, 3)
            centers = np.array([1, 2, 3])
            extra = {"full_list": model_find_parameters}
            return tracers, centers, extra

    finder = TestVFinder()

    # no copy
    voids = finder.find(box, box_copy=False)

    assert voids.tracers_in_voids_ == (1, 2, 3)
    np.testing.assert_array_equal(voids.centers_, [1, 2, 3])
    np.testing.assert_array_equal(voids.e_.full_list, [1, 2])
    assert voids.box is box

    # box copy
    voids = finder.find(box, box_copy=True)

    assert voids.tracers_in_voids_ == (1, 2, 3)
    np.testing.assert_array_equal(voids.centers_, [1, 2, 3])
    np.testing.assert_array_equal(voids.e_.full_list, [1, 2])
    assert voids.box is not box
    assert voids.box == box

