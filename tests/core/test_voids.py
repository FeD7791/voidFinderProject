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

"""Tests for the voidfindertk.core.voids."""

# =============================================================================
# IMPORTS
# =============================================================================

import os
import tempfile

import numpy as np

import pytest


from voidfindertk.core import voids
from voidfindertk.datasets import spherical_cloud


# =============================================================================
# TESTS
# =============================================================================


# @pytest.mark.parametrize("cleaner_method", ["overlap", "cbl"])
@pytest.mark.parametrize("cleaner_method", ["overlap"])
@pytest.mark.parametrize("radius_method", ["density", "extra", "volume"])
def test_voids(
    cleaner_method,
    radius_method,
    build_box_with_eq_voids,
    find_bubble_neighbors,
):

    rad = 30.0
    cloud = spherical_cloud.build_cloud()
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, rad=rad
    )

    dist, idx = find_bubble_neighbors(
        box=box,
        cloud_with_voids=cloud_with_voids,
        centers=centers,
        rad=rad * np.ones(len(centers)),
    )

    vds = voids.Voids(
        method="VoidFinderMethod",
        box=box,
        tracers_in_voids_=tuple(idx),
        centers_=centers,
        extra_={
            "radius": rad * np.ones(len(centers)),
            "properties": {
                "VoidVol": (4 / 3) * np.pi * (rad * np.ones(len(centers))) ** 3
            },
        },
    )

    with tempfile.TemporaryDirectory() as tempdir:

        vds.find_radius_and_clean(
            cleaner_method=cleaner_method,
            radius_method=radius_method,
            temporal_dir_path=os.path.abspath(tempdir),
        )


def test_voids_create_new_instance(
    build_box_with_eq_voids, find_bubble_neighbors
):
    rad = 30.0
    cloud = spherical_cloud.build_cloud()
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, rad=rad
    )

    dist, idx = find_bubble_neighbors(
        box=box,
        cloud_with_voids=cloud_with_voids,
        centers=centers,
        rad=rad * np.ones(len(centers)),
    )

    vds = voids.Voids(
        method="VoidFinderMethod",
        box=box,
        tracers_in_voids_=tuple(idx),
        centers_=centers,
        extra_={"radius": rad * np.ones(len(centers))},
    )
    parameters = {
        "method": "VoidFinderMethod",
        "box": box,
        "tracers_in_voids_": tuple(idx),
        "centers_": centers,
        "extra_": {"radius": rad * np.ones(len(centers))},
    }
    vds._create_new_instance(parameters=parameters)


@pytest.mark.parametrize("max_radius_search", [10, 30, 60, 70])
@pytest.mark.parametrize("delta_rad", [1, 5, 10, 20])
@pytest.mark.parametrize("n_jobs", [1, 4, 8])
def test_void_galaxy_corr(
    max_radius_search, delta_rad, n_jobs, generic_void_builder
):
    vds = generic_void_builder(rad=30.0)
    vds.void_galaxy_corr(
        max_radius_search=max_radius_search, delta_rad=delta_rad, n_jobs=n_jobs
    )


@pytest.mark.parametrize("wrong_n_jobs", [-11, 12.1, 0])
def test_void_galaxy_corr_wrong_n_jobs(wrong_n_jobs, generic_void_builder):
    vds = generic_void_builder(rad=30.0)
    with pytest.raises(ValueError, match="Invalid n_jobs value"):
        vds.void_galaxy_corr(
            max_radius_search=40, delta_rad=5, n_jobs=wrong_n_jobs
        )


@pytest.mark.parametrize("wrong_rad_search", [-1, 0])
def test_void_galaxy_corr_wrong_max_rad(
    wrong_rad_search, generic_void_builder
):
    vds = generic_void_builder(rad=30.0)
    with pytest.raises(ValueError, match="Invalid max_rad value"):
        vds.void_galaxy_corr(
            max_radius_search=wrong_rad_search, delta_rad=5, n_jobs=2
        )


@pytest.mark.parametrize("wrong_delta", [-1, 0])
def test_void_galaxy_corr_wrong_delta_r(wrong_delta, generic_void_builder):
    vds = generic_void_builder(rad=30.0)
    with pytest.raises(ValueError, match="Invalid delta_r value"):
        vds.void_galaxy_corr(
            max_radius_search=40, delta_rad=wrong_delta, n_jobs=2
        )
