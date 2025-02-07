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

""" test for voidfindertk.plot_acc

"""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np

import pytest

from voidfindertk.core import plot_acc
from voidfindertk.core import voids
from voidfindertk.datasets import spherical_cloud


def test_hist2d(mkbox):
    box = mkbox(seed=42, size=1000)
    plotter = plot_acc.BoxPlotter(box)

    ax = plotter.hist2d("x", "y")

    # Assert that the plot was created and labels are set correctly
    assert ax.get_xlabel() == f"x ({box.x.unit})"
    assert ax.get_ylabel() == f"y ({box.y.unit})"


def test_invalid_kind(mkbox):
    box = mkbox(seed=42, size=1000)
    plotter = plot_acc.BoxPlotter(box)

    # Test that a ValueError is raised for an invalid kind
    with pytest.raises(ValueError):
        plotter(kind="_private_method")


def test_void_size_function(build_box_with_eq_voids, find_bubble_neighbors):
    rad = 30.0
    cloud = spherical_cloud.build_cloud(n_points=100**3)
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, rad=rad
    )
    tracers_dist, tracers_ind = find_bubble_neighbors(
        box=box, cloud_with_voids=cloud_with_voids, centers=centers, rad=rad
    )
    v = voids.Voids(
        method="generic_method",
        box=box,
        tracers_in_voids_=tracers_ind,
        centers_=centers,
        extra_={"extra": ["extra_elements"]},
    )
    plotter = plot_acc.VoidPlotter(v)

    ax = plotter.void_size_function(
        vsf_kws={"radius": rad * np.ones(len(centers))}
    )

    # Assert that the plot was created and labels are set correctly
    # assert ax.get_xlabel() == r"$log_{10}(R)$R in unit"
    assert ax.get_ylabel() == r"$\frac{1}{V} \frac{dN_v}{dlnR_v}$"


def test_invalid_kind_void(build_box_with_eq_voids, find_bubble_neighbors):
    rad = 30
    cloud = spherical_cloud.build_cloud(n_points=100**3)
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, rad=rad
    )
    tracers_dist, tracers_ind = find_bubble_neighbors(
        box=box, cloud_with_voids=cloud_with_voids, centers=centers, rad=rad
    )
    v = voids.Voids(
        method="generic_method",
        box=box,
        tracers_in_voids_=tracers_ind,
        centers_=centers,
        extra_={"extra": ["extra_elements"]},
    )
    import ipdb

    ipdb.set_trace()
    plotter = plot_acc.VoidPlotter(v)

    # Test that a ValueError is raised for an invalid kind
    with pytest.raises(ValueError):
        plotter(kind="_private_method")
