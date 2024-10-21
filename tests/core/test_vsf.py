#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023 - 2024, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
#                            Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.

import numpy as np

from scipy.spatial.distance import cdist

from voidfindertk.core import box, vsf


def test_kind0_effective_radius_fixed_radius(
    build_spherical_void, build_cloud
):
    delta = -0.8
    rad = 50

    # Generate 1000 centers
    # np.random.seed(56)
    xyz_ = np.arange(0.0, 1000.0, rad + 10)
    centers = np.column_stack((xyz_, xyz_, xyz_))

    # Build cloud of points
    cloud = build_cloud(n_points=1000000, lmax=1000)

    # Get cloud with spherical voids
    cloud_with_voids = build_spherical_void(
        delta=delta, centers=centers, radii=rad, cloud=cloud
    )
    d = cdist(cloud_with_voids, np.array([[0.0, 0.0, 0.0]]))
    d = np.sort(np.ravel(d))

    # Get estimate of interparticle separation to the 10nt nearest neighbor.
    nn = 5
    separation = np.max(d[nn:] - d[:-nn])

    # Density threshold
    threshold = (
        (1 + delta) * len(cloud_with_voids) / (np.max(cloud_with_voids) ** 3)
    )
    # Transform the cloud to a box
    # box 1
    x, y, z = np.hsplit(cloud_with_voids, 3)
    b = box.Box(
        **{
            "x": x,
            "y": y,
            "z": z,
            # I don't care about this parameters.
            # But i have to fill them.
            "vx": x,
            "vy": y,
            "vz": z,
            "m": z,
        }
    )

    # Effective radius calculations
    er = vsf.effective_radius(
        centers=centers, n_neighbors=200, box=b, delta=delta, n_cells=64
    )

    # er has 4 parameters, errors, radius, tracers, densitymap
    error_kind0 = np.where(er.errors == 0)  # Normal (No error)
    # =========================================================================
    # ERROR KIND 0 TEST
    # =========================================================================

    # Get density values for this case:
    d0 = er.densities[
        error_kind0
    ]  # d0 is an array of arrays of density values.

    # d0_above_below = [[a0,b0],[a1,b1],...]
    # Here ai = number of denisties above threshold.
    # bi = number of densities below threshold.
    d0_above_below = np.array(
        [
            (len(np.where(_ > threshold)[0]), len(np.where(_ < threshold)[0]))
            for _ in d0
        ]
    )

    # Get the founded radii of each void
    rad0 = er.radius[error_kind0]

    # Assert some values are above threshold density and some are below.
    # .T[0] : Number of above values
    # .T[1] : Number of below values
    assert np.all(np.greater(d0_above_below.T[0], np.zeros(len(d0))))
    assert np.all(np.greater(d0_above_below.T[1], np.zeros(len(d0))))

    # The method should capture void radii with some presition.
    assert np.all(np.abs(rad0 - rad) < separation)
