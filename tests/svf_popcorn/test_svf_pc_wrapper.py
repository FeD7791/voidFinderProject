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

"""Test for voidfindertk.svf_popcorn._svf_pc_wrapper."""

# =============================================================================
# IMPORTS
# =============================================================================

import os
import pathlib
from unittest import mock


import numpy as np

import sh


from voidfindertk.svf_popcorn import _svf_pc_wrapper

# =============================================================================
# TESTS
# =============================================================================


def test_config_file_maker():
    path_this_directory_file = pathlib.Path(os.path.abspath(__file__))
    parameters = {
        "trsfile": "<trsfile>",
        "filefmt": "<filefmt>",
        "num_file": "<num_file>",
        "sphfile": "<sphfile>",
        "popfile": "<popfile>",
        "auxfiles": "<auxfiles>",
        "rawpopfile": "<rawpopfile>",
        "pairsfile": "<pairsfile>",
        "boxsize": "<boxsize>",
        "densth": "<densth>",
        "minradius": "<minradius>",
        "maxradius": "<maxradius>",
        "massmin": "<massmin>",
        "eps": "<eps>",
        # Directory path to place the file
        "path": str(path_this_directory_file.parents[0] / "vars.conf"),
    }

    # Create config file
    _svf_pc_wrapper.config_file_maker(**parameters)

    # Read the config file
    with open(parameters["path"], "r") as f:
        data = f.readlines()

    readed_parameters = dict(
        [
            (
                line.strip().split("=", 1)[0].strip(),
                line.strip().split("=", 1)[1].strip(),
            )
            for line in data[1:-1]
        ]
    )

    # Delete file
    sh.rm(parameters["path"])

    # work on values a little to compare them.
    parameters.pop("path")
    params_upper = {key.upper(): value for key, value in parameters.items()}

    assert readed_parameters == params_upper


def test_popcorn_svf_input_data_builder(mkbox):

    box = mkbox(seed=42)

    # Mocking the DataFrame constructor to return a mock dataframe
    mock_df_instance = mock.MagicMock()

    with mock.patch(
        "pandas.DataFrame", return_value=mock_df_instance
    ) as mock_dataframe:
        # Mock file path
        mock_file_path = mock.MagicMock(spec=pathlib.Path)

        # Call the function to test
        _svf_pc_wrapper.popcorn_svf_input_data_builder(
            box=box, file_path=mock_file_path
        )

        # Expected data that would be passed to the DataFrame constructor
        expected_data = np.array(
            [
                np.ravel(box.m),
                np.ravel(box.arr_.x),
                np.ravel(box.arr_.y),
                np.ravel(box.arr_.z),
                np.ravel(box.arr_.vx),
                np.ravel(box.arr_.vy),
                np.ravel(box.arr_.vz),
            ]
        ).T

    # Verify that DataFrame was created once with the expected data
    mock_dataframe.assert_called_once()

    args, kwargs = mock_dataframe.call_args
    assert np.array_equal(args[0], expected_data)  # Check data equality

    # Verify that the to_csv method was called once with the correct arguments
    mock_df_instance.to_csv.assert_called_once_with(
        mock_file_path, sep=" ", index=False, header=False
    )


def test_spherical_popcorn_void_finder():

    parameters = {
        "bin_path": pathlib.Path("<bin_path>"),
        "conf_file_path": "<conf_file_path>",
        "work_dir_path": ".",
        "cores": 1,
    }

    parameters_2 = {
        "bin_path": pathlib.Path("<bin_path>"),
        "conf_file_path": "<conf_file_path>",
        "work_dir_path": ".",
        "cores": None,
    }

    with mock.patch("sh.Command") as shcommand:
        _svf_pc_wrapper.spherical_popcorn_void_finder(**parameters)

    with mock.patch("sh.Command") as shcommand_2:
        _svf_pc_wrapper.spherical_popcorn_void_finder(**parameters_2)

    svf_popcorn = shcommand("svf", search_paths=[parameters["bin_path"]])
    shcommand("mpirun").assert_called_once_with(
        "-np",
        str(parameters["cores"]),
        "--bind-to",
        "core",
        str(svf_popcorn),
        f"config={parameters['conf_file_path']}",
        _cwd=parameters["work_dir_path"],
        _out=shcommand("ANY").call_args[1]["_out"],
        _err_to_out=True,
    )

    shcommand_2(
        "svf", search_paths=[parameters["bin_path"]]
    ).assert_called_once_with(
        f"config={parameters['conf_file_path']}",
        _cwd=parameters["work_dir_path"],
        _out=shcommand_2("ANY").call_args[1]["_out"],
        _err_to_out=True,
    )
