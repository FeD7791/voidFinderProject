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

"""Module to perform center calculatios based on tracers."""

# =============================================================================
# IMPORTS
# =============================================================================

from functools import reduce

import numpy as np

from voidfindertk.utils.box_to_grid import get_grispy_grid_from_box

# =============================================================================
# FUNCTIONS
# =============================================================================


def _get_xyz_tracers_from_indx(*, box, tracers_idx):
    """Find the (x,y,z) coordinates of tracer based on an index and the Box.

    Parameters
    ----------
        box : Box
            Object with tracer properties.

        tracers_idx : list
            Iterable with the indexes of the tracers
    """
    centers = np.column_stack((box.arr_.x, box.arr_.y, box.arr_.z))
    return [centers[idx] for idx in tracers_idx]


def center_calculator(box, tracers_in_voids, n_neighbors, threshold=0.8):
    """Calculates the centers of voids based on the positions of tracers \
    within each void.

    The function first constructs a grid from the provided 3D box and extracts
    the tracer positions from the specified voids. It then calculates the
    distance to the n_neighbors closest points for each tracer, determines the
    density of the tracers, and computes the weighted center of each void
    based on the tracer positions and densities.

    Border conditions are handled to ensure that the calculated centers fall
    within the specified box.

    Parameters
    ----------
    box : object
        A 3D box object representing the spatial domain. This is used to
        create a grid and define boundaries for the calculations.

    tracers_in_voids : list of list of int
        A list where each element corresponds to a void and contains the
        indices of tracers that are part of that void.

    n_neighbors : int
        The number of nearest neighbors to consider when calculating distances
        and densities.

    threshold : float, optional, default=0.8
        A threshold used to handle border conditions by adjusting tracer
        positions that are near the edge of the box. If a tracer's position
        exceeds the threshold times the box size, it is adjusted.

    Returns
    -------
    np.ndarray
        An array of shape (n_voids, 3) containing the weighted centers
        (x, y, z) of each void, where each center is calculated as the
        weighted average of tracer positions within the void, with weights
        determined by the density of tracers.

    Notes
    -----
    - The function assumes that the tracers are indexed by their positions in
    the box.
    - The grid for nearest neighbors is built using the
    `get_grispy_grid_from_box` function and requires an appropriate grid
    resolution (N_cells=64).
    - The density calculation is based on the inverse of the volume of a
    sphere with a radius equal to the farthest neighbor distance for each
    tracer.
    - The method also accounts for edge effects, ensuring centers are
    correctly positioned within the box.
    """
    # Build grispy grid from Box
    grid = get_grispy_grid_from_box(box=box, N_cells=64)
    # Get (x,y,z) positions from tracers in each void.
    # We get a list of positions of tracers for each void
    tinv = _get_xyz_tracers_from_indx(box=box, tracers_idx=tracers_in_voids)

    # Merge the list into a single array
    all_tracers_joined = reduce(
        lambda x, y: np.concatenate((x, y), axis=0), tinv
    )

    # The following is to account for border conditions:
    atj = np.where(
        all_tracers_joined > threshold * box.size(),  # Condition
        all_tracers_joined - box.size(),  # If condition condition fullfilled
        all_tracers_joined,
    )  # Else

    # Calculate the distances to n_neighbors to each of the tracers in the
    # array
    distances_result, nn = grid.nearest_neighbors(centres=atj, n=n_neighbors)

    radius = np.array([d[-1] for d in distances_result])
    density = 1 / (n_neighbors / ((4 / 3) * np.pi * radius**3))
    # Get the density as:
    # density = [n_neighbors/vol_sphere(r)]^(-1) (Inverse density actually)
    # r = farthest distane to n_neighbor

    # Secial index to iterate
    index_list = [0]
    [index_list.append(index_list[-1] + len(arr)) for arr in tinv]

    centers_ = np.array(
        [
            # For each set of k tracers in a single void: center =
            # (Sum i=0,k (x[i],y[i],z[i])*density[i])/ Sum(i=0,k density[k])
            # , i-est tracer in void
            np.sum(
                (value.T * density[index_list[idx] : index_list[idx + 1]]).T,
                axis=0,
            )
            / np.sum(density[index_list[idx] : index_list[idx + 1]])
            for idx, value in enumerate(tinv)
        ]
    )

    # Final version will account for border conditions
    # return np.where(centers_ < box.min(), centers_+box.size(), centers_)
    return centers_
