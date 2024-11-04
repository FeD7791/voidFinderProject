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
"""Module for process data obtained from SVF PopCorn void finder."""
# =============================================================================
# IMPORTS
# =============================================================================
import grispy as gsp

import numpy as np

import pandas as pd


def get_void_properties(*, popcorn_output_file_path):
    """
    Read void properties from a file and return them as a pandas dataframe.

    Parameters
    ----------
    popcorn_output_file_path : str
        Path to the file containing void properties. The file should
        contain lines where each line represents the properties of
        a void, with values separated by whitespace.

    Returns
    -------
    df_properties : DataFrame
        A pandas DataFrame with the properties of the void.

    Notes
    -----
    Each line in the file should contain values for the following
    attributes: 'id', 'r', 'x', 'y', 'z', 'delta_r'. The values
    should be separated by whitespace and appear in this order.
    """
    df_properties = pd.read_csv(
        popcorn_output_file_path,
        names=["id", "r", "x", "y", "z", "density_contrast"],
        delim_whitespace=True,
    )
    # return tuple(all_void_properties)
    return df_properties


def get_tracers_in_voids(*, box, popcorn_output_file_path):
    """
    Finds particles inside each void using Grispy.

    Parameters
    ----------
    box : object
        Box object (see box)

    popcorn_output_file_path : str
        Path to the file containing void properties. The file should
        have columns: 'id', 'rad', 'x', 'y', 'z', 'delta'. Each row
        represents a void with its properties including radius and
        center coordinates.

    Returns
    -------
    np.ndarray
        An array of indices of particles (refered to the Box object) that are
        within each void. The shape of the array depends on the implementation
        of `grid.bubble_neighbors`.

    Notes
    -----
    The function reads the void properties from the specified file,
    including void centers and radii. It then creates a Grispy grid
    using the provided particle coordinates and finds particles that
    lie within the voids based on the radius of each void.

    The `grid.bubble_neighbors` method is used to determine which
    particles fall within the specified distance (radius) from each
    void center.
    """
    # First get x,y,z array elements and combine them into a np.array(x,y,z):
    x = box.arr_.x
    y = box.arr_.y
    z = box.arr_.z
    xyz = np.column_stack((x, y, z))

    # Get radius and centers
    df = pd.read_csv(
        popcorn_output_file_path,
        delim_whitespace=True,
        names=["id", "rad", "x", "y", "z", "density_contrast"],
    )
    void_xyz = df[["x", "y", "z"]].to_numpy()
    void_rad = df["rad"].to_numpy()

    # Build grispy grid
    grid = gsp.GriSPy(xyz)
    # Set periodicity
    periodic = {
        0: (box.min(), box.max()),
        1: (box.min(), box.max()),
        2: (box.min(), box.max()),
    }
    grid.set_periodicity(periodic, inplace=True)

    # Get tracers
    dist, ind = grid.bubble_neighbors(void_xyz, distance_upper_bound=void_rad)
    return ind
