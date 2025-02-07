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


from voidfindertk.core.vfinder_abc import VoidFinderABC

# =============================================================================
# TESTS
# =============================================================================


def test_VFinderABC(mkbox):
    box = mkbox(seed=42)

    class TestVFinder(VoidFinderABC):
        def preprocess(self, box):
            return box

        def model_find(self, preprocess_parameters):
            return preprocess_parameters

        def build_voids(self, model_find_parameters):
            tracers_in_voids = [1, 2, 3]
            centers = [[0, 0, 0]]
            extra = {"extra": "extra"}

            return tracers_in_voids, centers, extra

    finder = TestVFinder()

    # no copy
    voids = finder.find(box, box_copy=False)

    assert voids.method == "TestVFinder"

    assert voids.box == box
    assert voids.tracers_in_voids_ == tuple([1, 2, 3])
    assert np.array_equal(voids.centers_, np.array([[0, 0, 0]]))
    assert voids.extra_ == {"extra": "extra"}
