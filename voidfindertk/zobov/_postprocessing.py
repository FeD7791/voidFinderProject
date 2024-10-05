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
"""Contains functions to parse output files from the ZOBOV void finder."""
# =============================================================================
# IMPORTS
# =============================================================================
import ctypes

import numpy as np


# =============================================================================
# FUNCTIONS
# =============================================================================
def parse_zones_in_void_output(
    *, executable_path, input_file_path, output_file_path
):
    """
    Parse tracers in zones output using a C library.

    Parameters
    ----------
    executable_path : str
        Path to the C library executable.
    input_file_path : str
        Path to the input file with tracers in zones data.
    output_file_path : str
        Path where the parsed output will be saved.

    Notes
    -----
    Uses ctypes to call a C function for parsing tracers in zones output.
    """
    # Get library
    clibrary = ctypes.CDLL(str(executable_path), mode=ctypes.RTLD_GLOBAL)

    # Input argtypes
    clibrary.process_files.argtypes = [ctypes.c_char_p, ctypes.c_char_p]

    # Call function
    clibrary.process_files(
        str(input_file_path).encode(), str(output_file_path).encode()
    )


def parse_tracers_in_zones_output(
    *, executable_path, input_file_path, output_file_path
):
    """
    Parse tracers in zones output using a C library.

    Parameters
    ----------
    executable_path : str
        Path to the C library executable.
    input_file_path : str
        Path to the input file with tracers in zones data.
    output_file_path : str
        Path where the parsed output will be saved.

    Notes
    -----
    Uses ctypes to call a C function for parsing tracers in zones output.
    """
    # Get library
    clibrary = ctypes.CDLL(str(executable_path), mode=ctypes.RTLD_GLOBAL)

    # Input argtypes
    clibrary.get_tracers_in_zones.argtypes = [ctypes.c_char_p, ctypes.c_char_p]

    # Call function
    clibrary.get_tracers_in_zones(
        str(input_file_path).encode(), str(output_file_path).encode()
    )


def _get_tracers_in_zones(*, tracers_in_zones_path):
    """
    Reads a file containing tracers in zones information and returns the\
    tracers associated with each zone.

    Parameters
    ----------
    tracers_in_zones_path : str
        Path to the file that contains information about tracers in each zone.


    Returns
    -------
    list of numpy.ndarray
        A list where each element is a NumPy array containing the tracer
        indices (as integers) corresponding to each zone. The list is ordered
        by the ascending index of the "CoreParticle" for each zone.

    """
    with open(tracers_in_zones_path, "r") as f:  # Read Parsed file
        zones_tracers = f.readlines()

    # List that will hold for each entrance, an array of member tracers
    tracers_in_zones = []

    for i, zp in enumerate(zones_tracers):
        # Deal with the format of the tracers in zones file
        if zp.startswith(" Nparticles"):
            tracers = np.array(zones_tracers[i + 2].split(" ")[:-1], dtype=int)
            # tracer[0] = CoreParticle index
            # tracers_in_zones is sorted in ascending order of CoreParticle
            # values.
            tracers_in_zones.append(tracers)

    return tracers_in_zones


def _get_zones_in_void(zones_in_void_file_path):
    """
    Gets zones belonging to voids.

    Read the output file containing zones in each void and returns an array
    maping zones to the void they belong.

    Parameters
    ----------
    zones_in_void_file_path : str
        Path to the output file containing zones in each void

    Returns
    -------
    list of numpy.ndarray
        A list of numpy arrays where the first element of each array is an
        index. The following elements are the zones inside the void, with
        the void index being the same as the first element of the array.
    """
    with open(zones_in_void_file_path, "r") as f:
        zones = f.readlines()
    zones_in_void = [np.array(zone.split(), dtype=int) for zone in zones[2:]]
    return zones_in_void


def get_tracers_in_voids(
    *, properties_dataframe, tracers_in_zones_path, zones_in_void_path
):
    """
    Adds tracer information about indexes of tracers inside each void.

    Parameters
    ----------
    properties_dataframe : pandas.DataFrame
        A DataFrame containing properties of different voids or zones. This
        DataFrame is obtained from the zobov output txt file with void proper-
        ties.

    tracers_in_zones_path : str
        Path to the file that contains tracers in zones data.

    zones_in_void_path : str
        Path to the file that defines zones associated with each void.

    Returns
    -------
    tuple
        A tuple containing:
        - properties_dataframe : pandas.DataFrame
            The updated DataFrame with one additional column:
            - 'Tracers_in_void': A column containing the combined tracers for
            each void.
        - tinv : list of numpy.ndarray
            A list of arrays, where each array contains the indices of tracers
            associated with each void. The list is sorted by the ascending or-
            der of void indices.

    Notes
    -----
    Indexes of tracers goes from [0,N] in zobov, so there is a direct mapping
    with the box object.

    """
    # Get the tracers in zones from the parsed file
    # tracers_in_zones is a dict where each key is CoreParticle index value.
    # keys is sorted in ascending order of CoreParticle values.
    tracers_in_zones = _get_tracers_in_zones(
        tracers_in_zones_path=tracers_in_zones_path
    )
    # Get the zones in each void from the parsed file
    zones_in_void = _get_zones_in_void(
        zones_in_void_file_path=zones_in_void_path
    )

    # Sort by CoreParticle ascending
    df = properties_dataframe.sort_values(by=["CoreParticle"])
    df["tinz"] = tracers_in_zones
    tinv = []

    for zone in zones_in_void:
        # Get all the indexes of tracers that are in each zone , then ...
        # ... merge them into single array that contains all the indexes of...
        # ... tracers that are inside a void that is form by the combination
        # ... of these zones.
        df_cut = df[df["FileVoid#"].isin(list(zone))]
        indx_tracers_in_void = np.concatenate(np.array(df_cut["tinz"]))
        # Index of tracers in zobov goes from [0,N]
        # You can confirm this by loking at the min value index for tracers
        # asociated with Void# = 1, min value goes from 0

        # Hierarchy: The void that will have this tracers assigned is the...
        # ...void with lowest "Void#"
        tinv.append((np.min(df_cut["Void#"]), indx_tracers_in_void))

    # Sort array based on the first value of tuples.

    tinv.sort(key=lambda x: x[0])
    tinv = [p[1] for p in tinv]

    properties_dataframe["Tracers_in_void"] = tinv
    # to test : len of tracers in Void elements should be the same as the
    # correlated Void#Part
    return properties_dataframe, tinv
