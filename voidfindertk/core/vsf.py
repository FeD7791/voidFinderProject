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

"""Module for effective radius and Void Size Calculation."""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np


# =============================================================================
# VSF
# =============================================================================


def void_size_function(
    radius,
    delta,
    box,
    n_step1=2,
    n_step2=10,
    scale_1_num_samples=7,
    scale_2_num_samples=2,
):
    """Computes the Void Size Function (VSF) for a set of radii.

    This function calculates the Void Size Function by creating a histogram of
    the logarithm of radii values, with different binning strategies for two
    scaling ranges. It returns the void size function as an instance of the
    `VSF` class, including the logarithm of the radius, density of voids, and
    associated scaling information.

    Parameters
    ----------
    radius : numpy.ndarray
        Array of radii values for which the Void Size Function is to be
        computed.
    box : Box
        An object representing the simulation box, used for volume and unit
        calculations.
    delta : float
        Density contrast used to scale the radii.
    n_step1 : int, optional
        Step size for the first scaling range. Default is 2.
    n_step2 : int, optional
        Step size for the second scaling range. Default is 10.
    scale_1_num_samples : int, optional
        Number of bins for the first scaling range. Default is 7.
    scale_2_num_samples : int, optional
        Number of bins for the second scaling range. Default is 2.

    Returns
    -------
    log_of_radius : numpy.ndarray
        Logarithm of the radius values for which the void size function is
        computed, scaled by the unit of the box.
    counts : numpy.ndarray
        Density of voids as a function of the logarithmic radius values.
    delta : float
        The density contrast value used for scaling the radii.

    Notes
    -----
    The function first computes the volume of the simulation box and the mean
    density of tracers, then calculates a scaling factor based on the given
    density contrast and the tracer number density. It constructs histogram
    bins over logarithmic radius values, computes the histogram and densities,
    and returns the cleaned, non-zero density values along with their
    corresponding logarithmic radii.

    The output `log_of_radius` values are scaled according to the box unit,
    which is derived from the `box.x.unit` attribute.
    """
    # Volume of the simulation (Box volume)
    vol = box.size() ** 3

    # Mean density of tracers
    rhomed = len(box) / vol

    # Number of tracers
    n = np.concatenate(
        [np.arange(6, 11, n_step1), np.arange(12, round(max(radius)), n_step2)]
    )

    # Scaling calculation
    # Radius in function of rhomed, tracers n and delta: R(rhomed, n , delta)
    # Radius = scl
    scl = np.log10((3 / (4 * np.pi) * n / rhomed / (1 + delta)) ** (1 / 3))

    mxlg = np.log10(max(radius))
    mxlg_scl_diff = mxlg - max(scl)

    # Histogram bins calculation
    bins = np.sort(
        np.concatenate(
            [
                scl[:-1],
                np.linspace(
                    max(scl),
                    max(scl) + mxlg_scl_diff * 0.5,
                    scale_1_num_samples,
                ),
                np.linspace(
                    max(scl) + mxlg_scl_diff * 0.5 + mxlg_scl_diff * 0.1,
                    mxlg,
                    scale_2_num_samples,
                ),
            ]
        )
    )

    # Histogram calculation
    h, bin_edges = np.histogram(np.log10(radius), bins=bins)

    # Mid values of the histogram
    mids = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    # Density calculation
    density = h / np.diff(bins) / vol

    # Clean values
    # Remove zeros
    index = np.where(density > 0.0)[0]  # Non zero elements index

    log_of_radius = mids[index] * box.x.unit
    counts = density[index]

    return log_of_radius, counts, delta
