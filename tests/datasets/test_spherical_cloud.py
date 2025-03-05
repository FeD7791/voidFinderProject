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

"""Tests for voidfindertk.datasets.spherical_cloud."""

# =============================================================================
# IMPORTS
# =============================================================================

import grispy as gsp

import numpy as np

import pytest

from voidfindertk.datasets import spherical_cloud

# =============================================================================
# FUNCTIONS
# =============================================================================


@pytest.mark.parametrize(
    "lmin, lmax, n_points, expected_shape, expected_range",
    [
        (0, 10, 5, (5, 3), (0, 10)),  # Test case 1: Small range and few points
        (
            1,
            100,
            4,
            (4, 3),
            (1, 100),
        ),  # Test case 2: Moderate range and few points
        (-5, 5, 10, (10, 3), (-5, 5)),  # Test case 3: Symmetric range around 0
        (
            0,
            1,
            2,
            (2, 3),
            (0, 1),
        ),  # Test case 4: Very small range, only 2 points
        (
            -100,
            100,
            3,
            (3, 3),
            (-100, 100),
        ),  # Test case 5: Larger range, 3 points
        (
            0,
            1000,
            100,
            (100, 3),
            (0, 1000),
        ),  # Test case 6: Larger range, 100 points
        (
            0,
            10,
            1000000,
            (1000000, 3),
            (0, 10),
        ),  # Test case 7: Large number of points
    ],
)
def test_build_cloud(lmin, lmax, n_points, expected_shape, expected_range):
    # Generate the cloud
    cloud = spherical_cloud.build_cloud(
        lmin=lmin, lmax=lmax, n_points=n_points, seed=2
    )

    # Check if the shape of the cloud is as expected
    assert cloud.shape == expected_shape

    # Check if all points are within the expected range
    assert np.all(cloud >= expected_range[0]) and np.all(
        cloud <= expected_range[1]
    )


@pytest.mark.parametrize(
    "log_mass_min, log_mass_max, cloud_shape, expected_range",
    [
        (
            10,
            12,
            (5, 4),
            (10, 12),
        ),  # Test case 1: Small cloud with small mass range
        (
            1,
            5,
            (3, 4),
            (1, 5),
        ),  # Test case 2: Small cloud with larger mass range
        (
            7,
            15,
            (10, 4),
            (7, 15),
        ),  # Test case 3: Larger cloud with a broad mass range
        (
            0,
            1,
            (100, 4),
            (0, 1),
        ),  # Test case 4: Large cloud with a narrow mass range
        (
            -5,
            5,
            (10, 4),
            (-5, 5),
        ),  # Test case 5: Negative mass values with a small cloud
    ],
)
def test_add_mass_to_cloud(
    log_mass_min, log_mass_max, cloud_shape, expected_range
):
    # Create a random cloud of points
    cloud = np.random.rand(
        cloud_shape[0], 3
    )  # cloud_shape[0] rows, 3 columns (x, y, z)

    # Add mass to the cloud
    updated_cloud = spherical_cloud.add_mass_to_cloud(
        cloud=cloud,
        log_mass_min=log_mass_min,
        log_mass_max=log_mass_max,
        seed=2,
    )

    # Check if the updated cloud has the expected shape (cloud_shape + 1 for
    # the mass column)
    assert updated_cloud.shape == (cloud_shape[0], 4)

    # Check if all log_mass values are within the specified range
    log_mass_values = updated_cloud[
        :, 3
    ]  # The log mass values should be in the 4th column
    assert np.all(log_mass_values >= expected_range[0]) and np.all(
        log_mass_values <= expected_range[1]
    )


@pytest.mark.parametrize("delta, radii", [(-0.9, 30), (-0.8, 50), (-0.7, 70)])
def test_build_spherical_void(delta, radii):

    centers = np.array(
        [
            [0.0, 0.0, 0.0],
            [100.0, 100.0, 100.0],
            [200.0, 200.0, 200.0],
            [300.0, 300.0, 300.0],
            [400.0, 400.0, 400.0],
            [500.0, 500.0, 500.0],
            [600.0, 600.0, 600.0],
        ]
    )

    cloud = spherical_cloud.build_cloud()

    # Density calculations
    cloud_volume = round(np.max(cloud) - np.min(cloud)) ** 3
    cloud_density = len(cloud) / cloud_volume
    density = (1 + delta) * cloud_density

    cloud_with_voids = spherical_cloud.build_spherical_void(
        delta=delta, centers=centers, radii=radii, cloud=cloud
    )

    # Build Grispy grid using cloud_with_void
    grid = gsp.GriSPy(cloud_with_voids)
    # Add periodicity
    periodic = {
        0: (np.min(cloud_with_voids), np.max(cloud_with_voids)),
        1: (np.min(cloud_with_voids), np.max(cloud_with_voids)),
        2: (np.min(cloud_with_voids), np.max(cloud_with_voids)),
    }
    grid.set_periodicity(periodic, inplace=True)
    # Perform buble search of neighbors around center
    distances, tracers = grid.bubble_neighbors(centers, radii)
    # Calculate densities
    densities = np.array(list(map(len, tracers))) / (
        (4 / 3) * np.pi * radii**3
    )
    # Check densities are the desired density
    assert np.all(np.less_equal(densities, density * np.ones(densities.shape)))


def test_build_spherical_overdensity():
    lmax = 1000
    rad = 30
    delta = 2
    cloud = spherical_cloud.build_cloud(lmax=lmax)
    # Build centers
    xyz_ = np.arange(0.0, lmax, rad + 10)
    centers = np.column_stack((xyz_, xyz_, xyz_))
    # Build cloud with overdensities
    c_w_o = spherical_cloud.build_spherical_overdensity(
        delta=delta, centers=centers, radii=rad, cloud=cloud
    )
    # Build grid
    periodic = {0: (0, lmax), 1: (0, lmax), 2: (0, lmax)}
    grid = gsp.GriSPy(c_w_o)
    grid.set_periodicity(periodic, inplace=True)
    # Perform buble search
    dist_, indx = grid.bubble_neighbors(centers, distance_upper_bound=rad)
    # Get mean density of universe
    rho_med = c_w_o.shape[0] / lmax**3

    number_tracers_in_voids = np.array(list(map(len, indx)))
    # densities
    dens = number_tracers_in_voids / ((4 / 3) * np.pi * rad**3)
    # Ideal density
    id_dens = rho_med * (delta + 1) * np.ones(len(dens))
    assert np.allclose(dens, id_dens, rtol=0.5)
