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

""" test for voidfindertk.svf_popcorn._svf_pc_postprocessing

"""

# =============================================================================
# IMPORTS
# =============================================================================

from unittest import mock

import numpy as np

from voidfindertk.datasets import spherical_cloud
from voidfindertk.svf_popcorn import _svf_pc_postprocessing

# =============================================================================
# TESTS
# =============================================================================


def test_get_tracers_in_voids(build_box_with_eq_voids, find_bubble_neighbors):
    # Notes: get_tracers_in_voids searches the tracers using grispy,
    # concretely bubble neighbors.

    # This test just test that the middle logic of the method does not mess
    # with the final result

    # This is an example of dataset provided
    rad = 50

    # Here a cloud of tracers is built
    cloud = spherical_cloud.build_cloud(n_points=100**3)

    # Here voids are created, removing tracers and leaving underdense regions.
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, rad=rad
    )

    # Using grispy, we search for tracers.
    dist, ind = find_bubble_neighbors(
        box=box, cloud_with_voids=cloud_with_voids, centers=centers, rad=rad
    )

    parameters = {"box": box, "popcorn_output_file_path": "<path>"}

    # We use mock so we do not go through the step of creating a file and
    # reading it with pandas
    mock_df = mock.MagicMock()
    mock_df.__getitem__.return_value.to_numpy.side_effect = [
        centers,
        rad * np.ones(len(centers)),
    ]

    with mock.patch("pandas.read_csv", return_value=mock_df) as mock_read_csv:
        ind_method = _svf_pc_postprocessing.get_tracers_in_voids(**parameters)

    # Both lists of indexes should be equal.
    assert all([np.all(a == b) for a, b in zip(ind_method, ind)])

    # Make sure you called pandas read_csv method.
    mock_read_csv.assert_called_once_with(
        parameters["popcorn_output_file_path"],
        **{
            "delim_whitespace": True,
            "names": ["id", "rad", "x", "y", "z", "density_contrast"],
        },
    )


def test_get_void_properties():
    mock_pd_read_csv = mock.MagicMock()
    parameters = {"popcorn_output_file_path": "<file_path>"}
    with mock.patch(
        "pandas.read_csv", return_value=mock_pd_read_csv
    ) as mock_read_csv:
        _svf_pc_postprocessing.get_void_properties(**parameters)

    mock_read_csv.assert_called_once_with(
        parameters["popcorn_output_file_path"],
        **{
            "names": ["id", "r", "x", "y", "z", "density_contrast"],
            "delim_whitespace": True,
        },
    )
