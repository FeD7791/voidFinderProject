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

"""Void-Galaxy Correlation Function Module."""

# =============================================================================
# IMPORTS
# =============================================================================

from joblib import Parallel, delayed

import numpy as np


from ..utils import box_to_grid


# =============================================================================
# FUNCTIONS
# =============================================================================


def _single_vgcf(centers, r_in, r_out, box, **grispy_kwargs):
    """
    Calculate the number density within a spherical shell of radius `r_in` \
    and `r_out` for a given set of centers.

    Parameters
    ----------
    centers : array_like, shape (N, 3)
        An array of 3D coordinates representing the positions of tracer
        particles. Each row represents a tracer's position in the form
        (x, y, z).

    r_in : float
        The inner radius of the shell. Should be smaller than `r_out`.

    r_out : float
        The outer radius of the shell. Should be larger than `r_in`.

    box : object
        A box object with periodic boundary conditions containing tracers
        properties.

    **grispy_kwargs : keyword arguments
        Additional keyword arguments passed to the `get_grispy_grid_from_box`
        function for grid construction.

    Returns
    -------
    np.ndarray, shape (M,)
        The number density of tracers per unit volume for each shell.

    """
    # Check r_in < r_out
    if r_in >= r_out:
        raise ValueError("Bad relation between radii, should be r_in < r_out")

    # Build periodic grid from box
    grid = box_to_grid.get_grispy_grid_from_box(box=box, **grispy_kwargs)

    # Find tracers between the shells
    shell_dist, shell_ind = grid.shell_neighbors(
        centers, distance_lower_bound=r_in, distance_upper_bound=r_out
    )

    # Get number of tracers within each shell
    n_tracers = np.array(list(map(len, shell_ind)))

    # Build Volumes
    v_arr = (4 / 3) * np.pi * (r_out**3 - r_in**3) * np.ones(len(n_tracers))

    # Data - Data value
    dd = n_tracers / v_arr

    # Random - Random value
    rr = len(box) / (box.size()) ** 3

    return np.sum(dd) / (len(dd) * rr) - 1


def vgcf_statistic(centers, box, max_rad, delta_r, n_jobs=1):
    """Calculates the Void-Galaxy Correlation Function.

    The correlation function is obtained as: vgcf = DD/RR - 1. The calculation
    is performed in concentric rings centered on void pre calculated centers.

    Parameters
    ----------
    centers : array_like, shape (N, 3)
        (x,y,z) coordinates of void centers.

    box : object
        A box object with periodic boundary conditions containing tracers
        properties.

    rad : array_like, shape (N,1)
        Array or radius values of the concentric rings.

    delta_r : float
        Width of the concentric rings

    n_jobs : int
        joblib parameter. The maximum number of concurrently running jobs.

    Returns
    -------
    vgcf : list
        Values of the vgcf.
    """
    # Check reasonable values
    if max_rad <= 0:
        raise ValueError("Invalid max_rad value")
    elif delta_r <= 0:
        raise ValueError("Invalid delta_r value")
    elif n_jobs < -1 or n_jobs == 0 or not isinstance(n_jobs, int):
        raise ValueError("Invalid n_jobs value")
    else:

        rad = np.arange(0, max_rad, delta_r)
        r_min = rad[:-1]
        r_max = rad[1:]

        # Set parallel mapping.
        parallel = Parallel(n_jobs=n_jobs, return_as="generator")

        vgcf = parallel(
            delayed(_single_vgcf)(  # Function
                centers=centers,  # Parameters
                box=box,
                r_in=rmin,  # loop variables
                r_out=rmax,
            )
            for rmax, rmin in zip(r_max, r_min)
        )

        return list(vgcf)
