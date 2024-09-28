#!/usr/bin/env python3
# =============================================================================
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
# Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
# =============================================================================

import grispy as gsp

import numpy as np

from voidfindertk.core import box, vsf


def test_effective_radius(mkbox):

    def build_cloud(*, seed=2, lmin=0, lmax=1000, n_points=100**3):
        """
        Builds a cloud of n_points , where: lmin < x,y,z < lmax.
        """
        np.random.seed(seed)
        # Create point cloud
        cloud = np.random.uniform(lmin, lmax, size=(n_points, 3))
        return cloud

    def build_spherical_void(
        delta,
        centers: np.ndarray,
        radii: float,
        cloud,
    ):
        """
        Build a cloud with voids with some density contrast delta at the
        centers provided by the user.
        """
        cloud_volume = round(np.max(cloud) - np.min(cloud)) ** 3
        cloud_density = len(cloud) / cloud_volume
        density_voids = (1 + delta) * cloud_density

        # Find neares neighbors
        grid = gsp.GriSPy(cloud)
        dist, index = grid.bubble_neighbors(
            centers, distance_upper_bound=radii
        )

        # Calculate right number of tracers so the void gets at the desired
        # density
        n = round(density_voids * (4 / 3) * np.pi * (radii**3))

        new_index = [
            i[: len(i) - n + 1] for i in index
        ]  # +1 to ensure we remove more particles so density of the void <
        # dens_voids

        # Remove from cloud the indicated indexes
        mask = np.ones(len(cloud), dtype=bool)
        for i in new_index:
            mask[i] = False

        cloud_with_voids = cloud[mask]
        return cloud_with_voids

    centers = np.array(
        [
            np.array([110, 110, 110]),
            np.array([120, 120, 120]),
            np.array([130, 130, 130]),
        ]
    )
    # Build cloud of points
    cloud = build_cloud(n_points=100000, lmax=180)
    # Get cloud with spherical voids
    rad = 30
    kloud = build_spherical_void(
        delta=-0.99, centers=centers, radii=rad, cloud=cloud
    )
    # Transform the cloud to a box
    x, y, z = np.hsplit(kloud, 3)
    b = box.Box(
        **{
            "x": x,
            "y": y,
            "z": z,
            # I don't care about this parameters
            "vx": x,
            "vy": y,
            "vz": z,
            "m": z,
        }
    )
    radii, tracers, density = vsf.effective_radius(
        centers=centers, n_neighbors=400, box=b, delta=-0.9
    )
    assert all(
        [
            abs(radii[0] - rad) < 0.2,
            abs(radii[1] - rad) < 0.2,
            abs(radii[2] - rad) < 0.2,
        ]
    )
