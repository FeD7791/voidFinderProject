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

from scipy.stats import linregress


from voidfindertk.core import center_finder
from voidfindertk.datasets.spherical_cloud import build_cloud

# =============================================================================
# TESTS
# =============================================================================


def test_center_calculator(build_box_with_eq_voids, find_bubble_neighbors):
    """
    Test to measure statistically that the right centers are found for this
    controlled example.
    """
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
    new_centers_ = center_finder.center_calculator(
        box=box,
        tracers_in_voids=indx,
        n_neighbors=10,
        threshold=0.9,
        n_jobs=1,
        batch_size=10,
    )

    # sort unordered centers
    new_centers_ = new_centers_[np.argsort(new_centers_.T[0])]
    # In this section, we remove the conditions that tracers are within box
    # This way comparison between cener locations are easier
    centers_for_test = []
    for center in new_centers_:
        counts = np.histogram(center, range=[0, box.size()])[0]
        if (counts[0] > 0) & (counts[-1] > 0):
            tracers_ = np.where(
                center > 0.8 * box.size(),  # Condition
                center - box.size(),  # If condition condition fullfilled
                center,  # Else
            )
            centers_for_test.append(tracers_)
        else:
            centers_for_test.append(center)

    centers_for_test = np.array(centers_for_test)

    # If for each coordinate center lie within a line with slope close to one
    # then center coordinates do agree. To check that we use the r coeficient.
    result_x = linregress(centers0.T[0], centers_for_test.T[0])
    result_y = linregress(centers0.T[1], centers_for_test.T[1])
    result_z = linregress(centers0.T[2], centers_for_test.T[2])

    assert (
        (result_x.rvalue > 0.9)
        & (result_y.rvalue > 0.9)
        & (result_z.rvalue > 0.9)
    )
