#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023 - 2024, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
#                            Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.

"""Module for performing void cleaning."""

# =============================================================================
# IMPORTS
# =============================================================================


import ctypes
import os
import pathlib
import shutil
import tempfile


import numpy as np

import pandas as pd

from scipy.spatial.distance import cdist


# =============================================================================
# CLEANER INTERFACE
# =============================================================================


def get_cleaner(*, cleaner_method):
    """Returns the appropriate cleaner function based on the specified method.

    Parameters
    ----------
    cleaner_method : str
        The cleaning method to use. Must be one of:
        - `"overlap"`: Uses `_overlap_cleaner`.
        - `"cbl"`: Uses `_cbl_cleaner`.

    Returns
    -------
    function
        The corresponding cleaner function.

    Raises
    ------
    ValueError
        If an invalid `cleaner_method` is provided.
    """
    if cleaner_method == "overlap":
        return _overlap_cleaner
    if cleaner_method == "cbl":
        return _cbl_cleaner
    else:
        ValueError(f"{cleaner_method} in not a valid cleaner method")


# =============================================================================
# DEFAULT CLEANER
# =============================================================================


def _overlap_cleaner(center, radius, **kwargs):
    """Removes overlapping voids based on their centers and radii by \
    comparing the distances between void centers and their combined radii.

    This function identifies pairs of voids that overlap (based on the sum of
    their radii) and removes the void with the smaller radius from each
    overlapping pair. If there are multiple overlaps, only the largest void is
    retained in each pair.

    Parameters
    ----------
    center : array_like, shape (n, 3)
        A 2D array containing the (x, y, z) coordinates of the centers of the
        voids. `n` represents the number of void centers.

    radius : array_like, shape (n,)
        A 1D array containing the radii of the voids. Each entry corresponds
        to the radius of the void with the same index in the `center` array.

    kwargs : keyword arguments, optional
        Additional parameters that can be passed for further customization.
        These are not currently used within the function, but they allow for
        future extensions.

    Returns
    -------
    tuple of (ndarray, ndarray)
        - center : ndarray, shape (m, 3)
            The coordinates of the void centers that do not overlap, where `m`
            is the number of non-overlapping voids.

        - radius : ndarray, shape (m,)
            The radii of the voids that do not overlap, corresponding to the
            void centers returned.

    Raises
    ------
    Warning
        If all voids overlap and have the same radius, a warning is raised
        stating that all voids overlap and none are retained.

    Notes
    -----
    This function computes a distance matrix between all pairs of void centers
    and checks for overlaps by comparing the distance between the centers with
    the sum of their radii. Voids with smaller radii in overlapping pairs are
    removed from the results.

    """
    size = len(center)

    # This matrix has zeros at diagonals
    dist = cdist(center, center)

    # Matrix with entrances r[i]+r[j]
    matrix0 = np.vstack([radius] * len(center))
    # np.eye(size)*2*radius will make diagonal elements zero in mat r[i]+r[j]
    matrix_ = matrix0 + matrix0.T - np.eye(size) * (2 * radius) - dist
    v_idx = np.unique(
        [
            np.where(col >= 0)[0][np.argmax(radius[np.where(col >= 0)[0]])]
            for col in matrix_
        ]
    )

    return center[v_idx], radius[v_idx]


# =============================================================================
# CBL CLEANER
# =============================================================================


def _cbl_cleaner(
    center,
    radius,
    box,
    *,
    # temp directory
    temporal_dir_path=".",
    clean_directory=True,
    # cbl parameters
    ratio=1.5,
    initial_radius=True,
    delta_r_min=10.0,
    delta_r_max=100.0,
    threshold=0.2,
    ol_crit="density_contrast",
    rescale=True,
    checkoverlap=True,
    **kwargs,
):
    """Sets working directories and input files to perform the CBL cleaning \
    using the function _cbl_cleaner_interface.

    Parameters
    ----------
    center : array_like
        An array of shape (n, 3), where n is the number of void centers,
        containing the (x, y, z) positions of the void centers.

    radius : array_like
        An array of shape (n,), where n is the number of void centers,
        containing the radius of the voids.

    box : object
        An object representing the box that contains properties of tracers.
        The structure of this object is assumed to be specific to the
        implementation and is passed to other functions for saving tracers and
        centers.

    temporal_dir_path : str, optional
        The path to the temporary directory where files will be stored during
        the process. Default is the current directory (`'.'`).

    clean_directory : bool, optional
        If True, the temporary directory and files will be removed after the
        function completes. Default is True.

    ratio : float, optional
        The ratio parameter passed to the CBL cleaning interface. Controls how
        the cleaning is performed. Default is 1.5.

    initial_radius : bool, optional
        If True, the initial radius will be considered during the cleaning
        process. Default is True.

    delta_r_min : float, optional
        Minimum delta radius for void cleaning.
        Default is 10.0.

    delta_r_max : float, optional
        Maximum delta radius for void cleaning.
        Default is 100.0.

    threshold : float, optional
        A threshold parameter controlling the cleaning procedure.
        Default is 0.2.

    ol_crit : str
        Criterion for the overlap step. Valid values are:
        [density_contrast, central_density]
        Default is density_contrast.

    rescale : bool, optional
        If True, rescaling is applied during the cleaning process.
        Default is True.

    checkoverlap : bool, optional
        If True, overlap checks will be performed during the cleaning.
        Default is True.

    kwargs : keyword arguments, optional
        Additional parameters passed to other functions in the workflow,
        typically to control additional options in the cleaning process. The
        exact parameters depend on the implementation of the specific cleaning
        functions.

    Returns
    -------
    tuple of (ndarray, ndarray)
        - centers : ndarray
            The cleaned (x, y, z) positions of the void centers, sorted by the
            norm of their positions.
        - radius : ndarray
            The cleaned radius values of the voids, sorted in the same order
            as the centers.

    Notes
    -----
    This function creates a temporary working directory, where it saves input
    files and calls the CBL cleaning interface to process the void centers and
    radii. After the cleaning is completed, the function reads the cleaned
    catalogue, removes the temporary directory (if requested), and returns the
    cleaned and sorted void centers and radii.

    """
    # Temporal Directory
    directory_path = pathlib.Path(
        tempfile.mkdtemp(suffix="cbl_cleaner", dir=temporal_dir_path)
    )
    cleaned_catalogue_path = (
        pathlib.Path(directory_path) / "cleaned_catalogue.txt"
    )
    # Setting variables
    input_tracers_path = directory_path / "input_tracers.txt"
    input_centers_path = directory_path / "input_xyz_rad.txt"

    # Overlap criterion
    if ol_crit == "density_contrast":
        _crit = True
    elif ol_crit == "central_density":
        _crit = False
    else:
        raise ValueError("Not a valid overlap Criterion variable")

    _save_xyz_tracers(box=box, path=input_tracers_path)
    _save_r_eff_center(centers=center, r_eff=radius, path=input_centers_path)
    _cbl_cleaner_interface(
        file_voids=input_centers_path,
        file_tracers=input_tracers_path,
        ratio=ratio,
        initial_radius=initial_radius,
        delta_r_min=delta_r_min,
        delta_r_max=delta_r_max,
        threshold=threshold,
        output_path=cleaned_catalogue_path,
        ol_crit=_crit,
        rescale=rescale,
        checkoverlap=checkoverlap,
    )
    centers, radius = _read_cleaned_catalogue(
        cleaned_catalogue_path=cleaned_catalogue_path
    )
    centers_norm = np.linalg.norm(centers, axis=1)
    # Remove temporal dir and files?
    if clean_directory:
        shutil.rmtree(directory_path)
    # Return arrays sorted ascending
    return centers[np.argsort(centers_norm)], radius[np.argsort(centers_norm)]


def _cbl_cleaner_interface(
    *,
    file_voids,
    file_tracers,
    ratio,
    initial_radius,
    delta_r_min,
    delta_r_max,
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
    the first 4 columns are going to be considered as inputs for
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
    to a selected range delta_r = delta_r = r_max - r_min for [r_min,r_max]
    are erased.

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
    file_voids_bytes = str(file_voids).encode("utf-8")
    file_tracers_bytes = str(file_tracers).encode("utf-8")

    delta_r_min_pointer = (ctypes.c_double * 1)(delta_r_min)
    delta_r_max_pointer = (ctypes.c_double * 1)(delta_r_max)
    output_path_bytes = str(output_path).encode("utf-8")

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
        ctypes.POINTER(ctypes.c_double),  # delta_r_min
        ctypes.POINTER(ctypes.c_double),  # delta_r_max
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
        delta_r_min_pointer,
        delta_r_max_pointer,
        threshold,
        output_path_bytes,
        ol_crit,
        rescale,
        checkoverlap,
    )


def _save_r_eff_center(*, centers, r_eff, path):
    """
    Save centers and effective radii to a tab-separated values file.

    Parameters
    ----------
    centers : array-like, shape (n_samples, 3)
        An array or list of coordinates representing the centers. Each entry
        should be a sequence of three values, corresponding to x, y, and z
        coordinates.

    r_eff : array-like, shape (n_samples,)
        An array or list of effective radii, where each value corresponds to
        the effective radius for the respective center.

    path : str
        The file path where the tab-separated values file will be saved.

    Returns
    -------
    None
        This function does not return any value.

    Notes
    -----
    The output file will have columns for x, y, z coordinates and r_eff,
    separated by tabs. The file will not include an index or header row.
    """
    df = pd.DataFrame(centers)
    df.columns = ["x", "y", "z"]
    df["r_eff"] = r_eff
    df.to_csv(path, index=False, header=False, sep="\t")


def _save_xyz_tracers(*, box, path):
    """
    Save x,y,z coordinates from box into a file.

    Save x, y, and z coordinates from a box object to a tab-separated values
    file.

    Parameters
    ----------
    box : object
        An object with attributes `x`, `y`, and `z`, each of which should have
        a `value` attribute that is an array-like sequence of coordinates.

    path : str
        The file path where the tab-separated values file will be saved.

    Returns
    -------
    None
        This function does not return any value.

    Notes
    -----
    The output file will have columns for x, y, and z coordinates, separated
    by tabs. The file will not include an index or header row.
    """
    x = box.x.value
    y = box.y.value
    z = box.z.value
    xyz = np.column_stack((x, y, z))
    df = pd.DataFrame(xyz)
    df.to_csv(path, index=False, header=False, sep="\t")


def _read_cleaned_catalogue(*, cleaned_catalogue_path):
    """Reads a void cleaned catalogue located in cleaned_catalogue_path and \
    returns the asociated center and radius.

    Parameters
    ----------
        cleaned_catalogue_path : str
            Path to te location of the catalgue file.

    Returns
    -------
        center, radius : tuple
        center : array of (x,y,z) positions of centers of voids
        radius : array of radius of voids.

    Notes
    -----
        param value in names is one of : [DensityContrast, CentralDensity]
        and this depends on ol_crit selected in parameters
        _cbl_cleaner_interface

    """
    df = pd.read_csv(
        cleaned_catalogue_path,
        delim_whitespace=True,
        names=["x", "y", "z", "rad", "param"],
    )
    center = np.array(df[["x", "y", "z"]])
    radius = np.array(df["rad"])
    return center, radius
