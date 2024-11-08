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

    # This is an example of dataset provided
    rad = 50
    cloud = spherical_cloud.build_cloud(n_points=100**2)
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, rad=rad
    )
    dist, ind = find_bubble_neighbors(
        box=box, cloud_with_voids=cloud_with_voids, centers=centers, rad=rad
    )

    parameters = {"box": box, "popcorn_output_file_path": "<path>"}

    mock_df = mock.MagicMock()
    mock_df.__getitem__.return_value.to_numpy.side_effect = [
        centers, rad * np.ones(len(box.arr_.x))
        ]

    with mock.patch("pandas.read_csv", return_value=mock_df) as mock_read_csv:
        ind_method = _svf_pc_postprocessing.get_tracers_in_voids(**parameters)
    import ipdb; ipdb.set_trace()
    assert np.all(ind_method == ind)
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
    with mock.patch("pandas.read_csv", return_value=mock_pd_read_csv):
        _svf_pc_postprocessing.get_void_properties(**parameters)

    mock_pd_read_csv.assert_called_once_with(
        parameters["popcorn_output_file_path"],
        **{
            "names": ["id", "r", "x", "y", "z", "density_contrast"],
            "delim_whitespace": True,
        },
    )
