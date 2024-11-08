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

""" tests for voidfindertk.datasets.spherical_cloud

"""

# =============================================================================
# IMPORTS
# =============================================================================

import grispy as gsp

import numpy as np

from voidfindertk.datasets import spherical_cloud

# =============================================================================
# FUNCTIONS
# =============================================================================


def test_cloud_with_overdensities():
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
