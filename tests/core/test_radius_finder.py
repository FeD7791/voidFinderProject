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

"""Test for voidfindertk.core.radius_finder."""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np

import pytest


from voidfindertk.core import radius_finder
from voidfindertk.datasets import spherical_cloud


# =============================================================================
# TESTS
# =============================================================================


@pytest.mark.parametrize(
    "delta, rad", list(zip(np.arange(-0.9, -0.7, 0.1), np.arange(20, 60, 10)))
)
def test_kind0_effective_radius(
    delta,
    rad,
    build_box_with_eq_voids,
    # get_first_neighbor_min_distance,
):

    # Beware: delta and n_neighbors are dependent variables. So if you increase
    # delta (this means to reduce the density), then you will have to increase
    # the n_neighbors. The dependence of each other is not calculated.

    # In this test, different values of delta are used, so that n_neighbors=200
    # are enough to find the right radius.

    # We want to check that radius are close enough to the real radius in an
    # spherical case.

    cloud = spherical_cloud.build_cloud()
    box, rho_threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        rad=rad, delta=delta, cloud=cloud
    )

    # Effective radius calculations
    eff_rad_obj = radius_finder.spherical_density_mapping(
        centers=centers, n_neighbors=400, box=box, delta=delta, n_cells=64
    )

    # er has 4 parameters, errors, radius, tracers, densitymap
    error_kind0 = np.where(eff_rad_obj.errors == 0)[0]  # Normal (No error)

    # Get density values for this case:
    d0 = eff_rad_obj.densities[
        error_kind0
    ]  # d0 is an array of arrays of density values.

    # d0_above_below = [[a0,b0],[a1,b1],...]
    # Here ai = number of denisties above rho_threshold.
    # bi = number of densities below rho_threshold.
    d0_above_below = np.array(
        [
            (
                len(np.where(_ > rho_threshold)[0]),
                len(np.where(_ < rho_threshold)[0]),
            )
            for _ in d0
        ]
    )
    # Get the founded radii of each void related to error kind 0.
    rad0 = eff_rad_obj.radius[error_kind0]
    # Finder should find ALL the voids for this kind of dataset.
    assert len(error_kind0) == len(eff_rad_obj.errors)

    # Assert some values are above rho_threshold density and some are below.
    # .T[0] : Number of above values
    # .T[1] : Number of below values
    assert np.all(d0_above_below.T[0] > np.zeros(len(d0)))
    assert np.all(d0_above_below.T[1] > np.zeros(len(d0)))

    # The method should capture void radii with some presition.
    # assert np.all(np.abs(rad0 - rad) < d_min)
    assert np.allclose(rad0, rad, atol=5)


@pytest.mark.parametrize("delta", np.arange(-0.94, -0.1, 0.05))
def test_kind3_effective_radius(delta, build_box_with_eq_voids):
    """Tests the output error UNDER_CRITICAL."""

    # Generate a cloud of tracers
    cloud = spherical_cloud.build_cloud()
    values = build_box_with_eq_voids(rad=50, cloud=cloud, delta=-0.95)
    centers = values[2]
    b = values[0]
    updated_cloud = values[3]

    with pytest.warns(
        RuntimeWarning,
    ) as record:
        # Effective radius calculations
        eff_rad_obj = radius_finder.spherical_density_mapping(
            centers=centers, n_neighbors=20, box=b, delta=delta, n_cells=64
        )

    errors_kind_3 = np.where(eff_rad_obj.errors == 3)[0]
    densities = eff_rad_obj.densities[errors_kind_3]
    max_dens = np.array(list(map(np.max, densities)))
    rho_threshold = 0.3 * len(updated_cloud) / (np.max(updated_cloud) ** 3)

    for message in record:
        assert "All values under critical Density for center " in str(message)

    # assert all densities last element is under rho_threshold
    assert np.max(max_dens) < rho_threshold


@pytest.mark.parametrize("delta", np.arange(-0.94, -0.1, 0.05))
def test_kind2_effective_radius(delta, build_box_with_eq_voids):
    """Tests the output error EXCEED_CRITICAL."""

    cloud = spherical_cloud.build_cloud()
    # Generate a cloud of tracers
    box, rho_threshold, centers, cloud = build_box_with_eq_voids(
        rad=20,
        cloud=cloud,
        delta=0.0,
    )
    # Perform Search
    eff_rad_obj = radius_finder.spherical_density_mapping(
        centers=centers, n_neighbors=50, box=box, delta=delta, n_cells=64
    )
    # Get all the kind 2 errors
    err_kind_2 = np.where(eff_rad_obj.errors == 2)[0]
    # Get the density map values regarding to those errors.
    dens = eff_rad_obj.densities[err_kind_2]
    # Get the min value of those errors.
    min_arr = np.array(list(map(np.min, dens)))
    # All densities should be above rho_threshold for this kind of error
    assert np.min(min_arr) > rho_threshold * (1 + delta)


def test_kind1_effective_radius(build_box_with_eq_voids, build_box_from_cloud):
    """Tests the output error MAYBE_NEAR_ANOTHER_VOID"""
    delta = -0.7
    # Generate a cloud of tracers
    cloud = spherical_cloud.build_cloud()
    b, t, centers, updated_cloud = build_box_with_eq_voids(
        rad=50, cloud=cloud, delta=-0.99
    )
    # Add overdensities
    cloud2 = spherical_cloud.build_spherical_overdensity(
        delta=2, centers=centers + 5, radii=10, cloud=updated_cloud
    )

    # Get box from cloud2
    b2 = build_box_from_cloud(cloud=cloud2)
    # Threshold
    rho_threshold = ((len(b2) / b2.size() ** 3)) * (1 + delta)

    # perform search
    eff_rad_obj = radius_finder.spherical_density_mapping(
        centers=centers, n_neighbors=100, box=b2, delta=delta, n_cells=64
    )

    err_kind_1 = np.where(eff_rad_obj.errors == 1)[0]
    dens = eff_rad_obj.densities[err_kind_1]

    down_threshold = np.array(
        [d[np.where(d < rho_threshold)[0]] for d in dens]
    )
    up_threshold = np.array([d[np.where(d > rho_threshold)[0]] for d in dens])
    last_element = np.array([d[-1] for d in dens])
    assert np.max(down_threshold) < rho_threshold
    assert np.min(up_threshold) > rho_threshold
    assert np.all(last_element < rho_threshold * np.ones(len(last_element)))


def test_effective_radius_class(build_box_with_eq_voids):

    delta = -0.9
    nn = 50
    n_cells = 64

    cloud = spherical_cloud.build_cloud()
    box, tracers, centers, updated_cloud = build_box_with_eq_voids(
        rad=30, cloud=cloud, delta=delta
    )

    effective_radius = radius_finder.spherical_density_mapping(
        centers=centers, n_neighbors=nn, box=box, delta=delta, n_cells=64
    )

    assert str(effective_radius) == (
        f"<effective_radius delta={delta} "
        f"n_neighbors={nn} n_cells={n_cells} | {len(centers)}/{len(centers)}>"
    )
    assert len(centers) == len(effective_radius)
    assert len(effective_radius.errors) == len(effective_radius.radius)
