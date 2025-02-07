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

"""
Test PopCorn
"""

# =============================================================================
# IMPORTS
# =============================================================================

from voidfindertk.datasets import spherical_cloud
from voidfindertk.popcorn import PopCorn

# =============================================================================
# TEST
# =============================================================================


def test_popcorn_working_example(build_box_with_eq_voids):

    delta = -0.9
    cloud = spherical_cloud.build_cloud(lmin=0, lmax=1000)
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, delta=delta, rad=30
    )
    model = PopCorn(
        workdir_clean=True,
        boxsize=1000,
        densth=delta,
        minradius=5,
        maxradius=100,
        shot_noise_threshold=15,
    )

    model.find(box=box)


def test_popcorn_box_tracers_mapping(build_box_with_eq_voids, find_tracers):

    delta = -0.8
    cloud = spherical_cloud.build_cloud(lmin=0, lmax=1000, n_points=100**3)
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, delta=delta, rad=30, log_mass_min=10, log_mass_max=50
    )

    model = PopCorn(
        workdir_clean=True,
        boxsize=1000,
        densth=delta,
        minradius=5,
        maxradius=100,
        shot_noise_threshold=15,
    )

    void = model.find(box=box)
    df = void.extra_.voids
    df["tracers2"] = void.tracers
    assert all(df.tracers == df.tracers2)
