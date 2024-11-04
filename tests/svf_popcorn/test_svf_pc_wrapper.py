#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023 - 2024, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
#                            Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.

# =============================================================================
# IMPORTS
# =============================================================================
import configparser
import os
import pathlib
from unittest import mock

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
        "path": str(path_this_directory_file / "vars.conf"),
    }
    file_parser = mock.MagicMock()
    with mock.patch("configparser.ConfigParser", return_value=file_parser):
        _svf_pc_wrapper.config_file_maker(**parameters)

    _svf_pc_wrapper.config_file_maker(**parameters)
    config = configparser.ConfigParser()
    config.read(parameters["path"])
    config_input_params = config["INPUT_PARAMS"]
    # Delete file
    sh.rm(parameters["path"])

    file_parser.assert_called_once_with(**parameters)
    assert all(list(parameters.values()) == list(config_input_params.values()))


def test_popcorn_svf_input_data_builder(mkbox):
    parameters = {"box": mkbox(seed=42), "file_path": "<filepath.txt>"}
    to_csv_parameters = {"sep": " ", "index": False, "header": False}

    mock_to_csv = mock.MagicMock()
    mock_dataframe = mock.MagicMock()
    with mock.patch("pandas.DataFrame", spec=True) as dataframe:
        dataframe.return_value = mock_dataframe
        mock_dataframe.to_csv.return_value = mock_to_csv
        _svf_pc_wrapper.popcorn_svf_input_data_builder(**parameters)

    dataframe.assert_called_once_with(**parameters)
    dataframe.to_csv.assert_called_once_with(
        parameters["file_path"], **to_csv_parameters
    )


def test_spherical_popcorn_void_finder():
    parameters = {
        "mpi_flags": None,
        "bin_path": "<bin_path>",
        "conf_file_path": "<conf_file_path>",
        "work_dir_path": "<work_dir_path>",
    }
    mock_sh = mock.MagicMock()
    mock_chdir = mock.MagicMock()
    with mock.patch("sh.Command", return_value=mock_sh):
        with mock.path(
            "voidfindertk.utils.context_managers.chdir",
            return_value=mock_chdir,
        ):
            _svf_pc_wrapper.spherical_popcorn_void_finder(**parameters)

    mock_sh.assert_called_once_with("svf", **{"bin_path": ["<bin_path>"]})
    mock_chdir.assert_called_once_with(parameters["work_dir_path"])
