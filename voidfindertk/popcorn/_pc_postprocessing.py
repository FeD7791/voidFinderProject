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

"""Module for process data obtained from PopCorn void finder."""

# =============================================================================
# IMPORTS
# =============================================================================


import numpy as np

import pandas as pd

# =============================================================================
# FUNCTIONS
# =============================================================================


def _to_float(line):
    line = line.strip().split()
    return list(map(np.float32, line))


def get_properties(filename):
    """Parse the output file from the PopCorn Void Finder and extract \
    properties of voids and spheres.

    This function reads the output file, processes its contents, and returns
    two pandas DataFrames: one containing the properties of voids, and the
    other containing the properties of spheres. The function assumes a specific
    format for the input file, where each line contains space-separated values
    corresponding to the properties of voids and spheres.

    Parameters
    ----------
    filename : str
        The path to the output file from the PopCorn Void Finder.

    Returns
    -------
    voids : pandas.DataFrame
        A DataFrame containing the properties of voids. The columns are:
        - 'id': The identifier of the void (int).
        - 'n_mem': The number of members in the void (int).
        - 'volume': The volume of the void (float).
        - 'n_part': The number of particles in the void (int).
        - 'flag': A flag indicating void properties (int).
        - 'tracers': A list of tracers (1D array).

    sphere : pandas.DataFrame
        A DataFrame containing the properties of spheres. The columns are:
        - 'x': The x-coordinate of the sphere's center (float).
        - 'y': The y-coordinate of the sphere's center (float).
        - 'z': The z-coordinate of the sphere's center (float).
        - 'radius': The radius of the sphere (float).
        - 'fracvol': The fractional volume of the sphere (float).
        - 'level': The level of the sphere (int).
        - 'id': The identifier of the void to which the sphere belongs (int).

    Example
    -------
    voids, sphere = get_properties('popcorn_output.txt')
    """
    with open(filename, "r") as f:
        data = f.readlines()  # Get all the lines in file

    # map each string in line to a list of spaced separeted values
    data_outputs = list(map(_to_float, data[1:]))

    voids = []
    spheres = []
    tracers = []
    idx_spheres = []

    for idx, out in enumerate(data_outputs):
        if len(out) == 5:
            voids.append(out)
            spheres.append(data_outputs[idx + 1 : idx + 1 + int(out[1])])
            tracers.append(
                np.ravel(
                    np.array(
                        data_outputs[
                            idx
                            + 1
                            + int(out[1]) : idx
                            + 1
                            + int(out[1])
                            + int(out[3])
                        ],
                        dtype=int,
                    )
                )
            )
            idx_spheres += [int(out[0])] * int(out[1])

    void_properties_names = ["id", "n_mem", "volume", "n_part", "flag"]
    sphere_properties_names = ["x", "y", "z", "radius", "fracvol", "level"]

    voids = pd.DataFrame(voids, columns=void_properties_names)
    voids = voids.astype({"id": "int", "n_mem": "int", "n_part": "int"})
    voids["tracers"] = tracers

    sphere = pd.DataFrame(np.vstack(spheres), columns=sphere_properties_names)
    sphere["id"] = idx_spheres

    return voids, sphere
