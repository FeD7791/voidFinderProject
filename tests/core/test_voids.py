#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023 - 2024, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
#                            Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.

# =============================================================================
# IMPORTS
# =============================================================================


import numpy as np

import pytest


from voidfindertk.core import voids
from voidfindertk.datasets import spherical_cloud


# =============================================================================
# TESTS
# =============================================================================
# @pytest.mark.parametrize(
#     "cleaner_method, radius_method",
#     list(zip(["overlap","cbl"],["default","extra","volume"])))
@pytest.mark.parametrize("cleaner_method", ["overlap", "cbl"])
@pytest.mark.parametrize("radius_method", ["default", "extra", "volume"])
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

    vds.find_radius_and_clean(
        cleaner_method=cleaner_method, radius_method=radius_method
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
