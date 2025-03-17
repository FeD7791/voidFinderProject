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


"""Tests for Void-Galaxy Correlation Function Module."""


# =============================================================================
# IMPORTS
# =============================================================================


import numpy as np

import pytest


from voidfindertk.core import vgcf
from voidfindertk.datasets import spherical_cloud


# =============================================================================
# TESTS
# =============================================================================


@pytest.mark.parametrize("delta", [-0.6, -0.7, -0.8, -0.9])
@pytest.mark.parametrize("delta_r", [2, 4, 6, 10, 20])
def test_vgcf(delta, delta_r, build_box_with_eq_voids):

    cloud = spherical_cloud.build_cloud()
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        delta=delta, rad=30.0, cloud=cloud, sep=100
    )
    v_g_cf = vgcf.vgcf_statistic(
        centers, box, max_rad=80, delta_r=delta_r, n_jobs=8
    )
    v_g_cf = np.array(v_g_cf)

    vgcf_abs = np.abs(v_g_cf)

    nn = int(np.rint(8 * len(v_g_cf) / 10))

    assert np.all(vgcf_abs[nn:] < 0.2)


@pytest.mark.parametrize("wrong_n_jobs", [-11, 12.1, 0])
def test_void_galaxy_corr_wrong_n_jobs(wrong_n_jobs, build_box_with_eq_voids):
    cloud = spherical_cloud.build_cloud()
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud
    )

    with pytest.raises(ValueError, match="Invalid n_jobs value"):
        vgcf.vgcf_statistic(
            centers, box, max_rad=40, delta_r=5, n_jobs=wrong_n_jobs
        )


@pytest.mark.parametrize("wrong_rad_search", [-1, 0])
def test_void_galaxy_corr_wrong_max_rad(
    wrong_rad_search, build_box_with_eq_voids
):
    cloud = spherical_cloud.build_cloud()
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud
    )

    with pytest.raises(ValueError, match="Invalid max_rad value"):
        vgcf.vgcf_statistic(
            centers, box, max_rad=wrong_rad_search, delta_r=5, n_jobs=4
        )


@pytest.mark.parametrize("wrong_delta", [-1, 0])
def test_void_galaxy_corr_wrong_delta_r(wrong_delta, build_box_with_eq_voids):
    cloud = spherical_cloud.build_cloud()
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud
    )

    with pytest.raises(ValueError, match="Invalid delta_r value"):
        vgcf.vgcf_statistic(
            centers, box, max_rad=40, delta_r=wrong_delta, n_jobs=4
        )
