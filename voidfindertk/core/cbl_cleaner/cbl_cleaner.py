#!/usr/bin/env python3
# =============================================================================
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
# Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
# =============================================================================
"""Wrapper of the Cosmo Bologna Lib clean_void_catalogue."""
import ctypes
import os
import pathlib

import numpy as np


def cbl_cleaner(
    *,
    file_voids,
    file_tracers,
    ratio,
    initial_radius,
    delta_r,
    threshold,
    output_path,
    ol_crit,
    rescale,
    checkoverlap,
):
    """
    Performs the CBL cleaning over a void catalogue.

    The cleaning is done by providing radius and center coordinates of a
    void catalogue provided by a void finder method.

    Parameters
    ----------
    file_voids : str
    Path to the file that holds a void catalogue. By default
    the first 4 columns are going tobe considered as inputs for
    x,y,z,r_eff where x,y,z are the barycentre of the void (see
    get_center_and_radii).

    file_tracers : str
    Path to the file that holds the input tracers. By default the first 3
    columns of the file are considered as inputs x,y,z refering to the
    positions of each tracer.

    ratio : float (0 < ratio < 1)
    Distance from the void centre at which the density contrast is
    evaluated in units of the void radius. Ex: ratio = 0.1
    =>  10% of the void radius lenght

    initial_radius : bool
    If true erase voids with effective radii outside a given range delta_r.

    delta_r : list
    Interval of accepted radii. Voids whose effective radius do not belong
    to a selected range delta_r = [r_min,r_max] are erased.

    threshold : float
    Erase voids with underdensities higher than a given threshold.


    output_path : str
    Path to the output cleaned catalogue.

    ol_crit : bool
    The criterion for the overlap step.
        True : The void with the lower density constrast is rejected.
        False : The void with the higher central density is rejected.

    Notes
    -----
    This function calculates the central density and the density contrast
    automatically using the ratio input variable.

    The central density (in units of the average density) is
    computed as the density of a sphere centred in the void centre and
    with radius R = ratio * R_eff.

    The density contrast is the ratio between the central density and the
    density within the sphere centred in the void centre and with radius:
    R = R_eff
    """
    # Prepare variables
    file_voids_bytes = file_voids.encode("utf-8")
    file_tracers_bytes = file_tracers.encode("utf-8")
    delta_r_array = np.array(delta_r, dtype=np.float64)
    delta_r_ctypes = delta_r_array.ctypes.data_as(
        ctypes.POINTER(ctypes.c_double)
    )
    output_path_bytes = output_path.encode("utf-8")

    # Path to the library
    path = os.path.dirname(os.path.abspath(__file__))
    clibrary = ctypes.CDLL(
        str(pathlib.Path(path) / "libcleaner.so"), mode=ctypes.RTLD_GLOBAL
    )

    # Input arguments
    clibrary.process_catalogues.argtypes = [
        ctypes.c_char_p,  # file_voids
        ctypes.c_char_p,  # file_tracers
        #  Parameters of Cleaner
        ctypes.c_double,  # ratio
        ctypes.c_bool,  # initial_radius
        ctypes.POINTER(ctypes.c_double),  # delta_r (pointer to double array)
        ctypes.c_int,  # delta_r_size
        ctypes.c_double,  # threshold
        ctypes.c_char_p,  # output_path
        ctypes.c_bool,  # ol_criterion
        ctypes.c_bool,  # rescale
        ctypes.c_bool,  # checkoverlap
    ]

    # Output Arguments
    clibrary.process_catalogues.restype = None

    # Call the C++ function
    clibrary.process_catalogues(
        file_voids_bytes,
        file_tracers_bytes,
        ratio,
        initial_radius,
        delta_r_ctypes,
        len(delta_r),
        threshold,
        output_path_bytes,
        ol_crit,
        rescale,
        checkoverlap,
    )
