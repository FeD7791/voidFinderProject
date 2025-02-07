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
import struct


import numpy as np


# =============================================================================
# FILE PARSERS
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


# =============================================================================
# FILE READERS
# =============================================================================


def _get_tracers_in_zones(*, tracers_in_zones_path):
    """
    Reads a file containing tracers in zones information and returns the \
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

        Void# in txt file is the same as first element of each array.
    """
    with open(zones_in_void_file_path, "r") as f:
        zones = f.readlines()
    zones_in_void = [np.array(zone.split(), dtype=int) for zone in zones[2:]]
    return zones_in_void


# =============================================================================
# MAIN POSTPROCESS
# =============================================================================


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


# =============================================================================
# CENTERS
# =============================================================================


def get_center_method(center_method):
    """Perform method selection based on string reference."""
    if center_method == "barycentre":
        return _centers_barycentre_method
    if center_method == "core_particle":
        return _center_core_particle_method


# =============================================================================
# CENTERS 1) By Core Particle
# =============================================================================


def _center_core_particle_method(*, properties_df, box):
    """
    Extracts the coordinates of core particles from a simulation box.

    This function retrieves the (x, y, z) coordinates from the given `box`
    and selects only the core particles based on indices provided in
    `properties_df["CoreParticle"]`.

    Parameters
    ----------
    properties_df : pandas.DataFrame
        A DataFrame containing at least a "CoreParticle" column, which
        holds indices indicating the core particles to extract.
    box : object
        An object with an attribute `arr_` containing `x`, `y`, and `z`
        attributes, representing spatial coordinates.

    Returns
    -------
    numpy.ndarray
        A `(N, 3)` array containing the (x, y, z) coordinates of the selected
        core particles.

    """
    x = box.arr_.x
    y = box.arr_.y
    z = box.arr_.z
    indx = properties_df["CoreParticle"]
    xyz = np.column_stack((x, y, z))
    return xyz[indx]


# =============================================================================
# CENTERS 2) By Baricentre
# =============================================================================


def _read_volume_file(*, filename):
    """
    Read data from the binary volume file and return it as a NumPy array.

    Parameters
    ----------
    filename : str
        The path to the binary volume file containing the volume data. This
        file contains the volumes of the voronoi cells where each particle is
        in.

    Returns
    -------
    volumes : numpy.ndarray
        A 1-D NumPy array of type `np.float32` containing the volume data read
        from the binary file.

    Notes
    -----
    The indexes of the returned array are directly related to the tracers
    index in the box object.

    Examples
    --------
    >>> volumes = read_volume_file(filename='voloutput_vozinit.dat')
    >>> print(volumes)
    [1.234 5.678 9.101]
    """
    with open(filename, "rb") as f:
        # Read number of tracers
        number_voids = struct.unpack("i", f.read(4))[0]
        # Read volumes
        volumes = np.zeros(number_voids, dtype=np.float32)
        for i in range(number_voids):
            volume = struct.unpack("d", f.read(8))[0]
            volumes[i] = np.float32(volume)
    return volumes


def _get_tracers_xyz(*, box):
    """
    Gets the x,y,z coordinates of the tracers inside the box.

    Parameters
    ----------
    box : Object Box
        Object that holds the properties of the input tracers

    Returns
    -------
    xyz_arr : numpy.array
        Array of N rows, 3 cols where each row is the x,y,z position of a
        tracer.
    """
    tracer_x = box.x.value
    tracer_y = box.y.value
    tracer_z = box.z.value
    xyz_arr = np.stack([tracer_x, tracer_y, tracer_z], axis=1)
    return xyz_arr


def _calculate_barycentre(*, tracers_xyz, tracers, tracer_volumes):
    """
    Calucates the barycentre of a single void.

    A void is composed of several voronoi cells. Each particle is inside a
    unique voronoi cell. The barycentre of a particular void iscalculated as
    the weigthed sum of the particles position times the volume of its voronoi
    cell divided the total volume of the void.

    Parameters
    ----------
    tracers_xyz : numpy.array
        Array of [x,y,z] tracers positions
    tracers : list
        List of indexes of each tracer that belongs to a particular void.
    tracer_volumes : numpy.array
        Array [v1,v2,v3...] of voronoi cell volumes of each particle.

    Returns
    -------
    center : numpy.array
        [X,Y,Z] coordinates of the center of the void.


    """
    tracer_volumes = tracer_volumes[tracers]
    arr = tracers_xyz[tracers]
    center = np.average(arr, weights=tracer_volumes, axis=0)
    return center


def _centers_barycentre_method(
    *, tracers_volumes_file_path, tracers_in_voids, box
):
    """
    Calculates the center by barycentre method of each void.

    Parameters
    ----------
    tracer_volumes : list
        List of the voronoi cell volumes holding each particle.
    tracers_in_voids : list
        List of particles within voids.
    box : Object
        Box Object that holds information about the tracers data.

    Returns
    -------
    centers : list
        List of (x,y,z) coordinates of each void center.

    """
    # Get volumes
    tracer_volumes = _read_volume_file(filename=tracers_volumes_file_path)
    # Get tracers xyz coords
    tracers_xyz = _get_tracers_xyz(box=box)

    # Get centers
    centers = []
    for tracers_in_void in tracers_in_voids:
        center = _calculate_barycentre(
            tracers_xyz=tracers_xyz,
            tracers=tracers_in_void,
            tracer_volumes=tracer_volumes,
        )
        centers.append(center)
    centers = np.array(centers)
    return centers
