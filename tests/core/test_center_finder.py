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

"""Test for center finder methods."""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np


from voidfindertk.core import center_finder
from voidfindertk.datasets.spherical_cloud import build_cloud


# =============================================================================
# TESTS
# =============================================================================


def test_center_calculator(build_box_with_eq_voids, find_bubble_neighbors):

    rad = 30
    cloud = build_cloud()
    box, threshold, centers0, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud
    )

    # First find the tracers within each void.
    dist, indx = find_bubble_neighbors(
        box=box, cloud_with_voids=cloud_with_voids, centers=centers0, rad=rad
    )

    # Find centers
    new_centers_1 = center_finder.center_calculator(
        box=box, tracers_in_voids=indx, n_neighbors=10, threshold=0.8
    )

    # This condition, is not suficient, you should have a precalculated
    # condition to check it.
    assert np.all(
        np.linalg.norm(new_centers_1[1:] - centers0[1:], axis=1) < 50
    )
