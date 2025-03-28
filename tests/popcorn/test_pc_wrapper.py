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

from unittest import mock

from voidfindertk.popcorn import _pc_wrapper

# =============================================================================
# TESTS
# =============================================================================


def test_popcorn_void_finder():
    parameters1 = {
        "bin_path": "bin_path",
        "conf_file_path": "conf_file_path",
        "work_dir_path": ".",
        "cores": 1,
    }

    parameters2 = {
        "bin_path": "bin_path",
        "conf_file_path": "conf_file_path",
        "work_dir_path": ".",
        "cores": None,
    }

    with mock.patch("sh.Command") as shcommand:
        _pc_wrapper.popcorn_void_finder(**parameters1)

    with mock.patch("sh.Command") as shcommand2:

        _pc_wrapper.popcorn_void_finder(**parameters2)

    popcorn = shcommand2("popcorn", search_paths=["bin_path"])

    shcommand("mpirun").assert_called_with(
        "-np",
        str(parameters1["cores"]),
        "--bind-to",
        "core",
        str(shcommand("popcorn", search_paths=[parameters1["bin_path"]])),
        f"config={parameters1['conf_file_path']}",
        _cwd=parameters1["work_dir_path"],
        _out=shcommand("mpirun").call_args[1]["_out"],
        _err_to_out=True,
    )

    popcorn.assert_called_with(
        f"config={parameters2['conf_file_path']}",
        _cwd=parameters2["work_dir_path"],
        _out=shcommand2("mpirun").call_args[1]["_out"],
        _err_to_out=True,
    )
    # Had to take out lambda function explicitly otherwise is imposible to get
    # the same reference


def test_compute_intersects():
    parameters = {
        "bin_path": "bin_path",
        "conf_file_path": "conf_file_path",
        "work_dir_path": ".",
    }

    with mock.patch("sh.Command") as shcommand:
        _pc_wrapper.compute_intersects(**parameters)

    shcommand(
        "compute_intersecs", search_paths=[parameters["bin_path"]]
    ).assert_called_once_with(
        f"config={parameters['conf_file_path']}",
        _cwd=parameters["work_dir_path"],
        _out=shcommand("ANY").call_args[1]["_out"],
        _err_to_out=True,
    )


def test_clean_duplicates():
    parameters = {
        "bin_path": "bin_path",
        "conf_file_path": "conf_file_path",
        "work_dir_path": ".",
    }

    with mock.patch("sh.Command") as shcommand:
        _pc_wrapper.clean_duplicates(**parameters)

    shcommand(
        "clean_duplicates", search_paths=[parameters["bin_path"]]
    ).assert_called_once_with(
        f"config={parameters['conf_file_path']}",
        _cwd=parameters["work_dir_path"],
        _out=shcommand("ANY").call_args[1]["_out"],
        _err_to_out=True,
    )


def test_read_and_modify_config():
    # mock_config = mock.create_autospec(configparser.ConfigParser)
    parameters = {
        "config_file_path": "config_file_path",
        "section": "section",
        "parameter": "parameter",
        "new_value": "new_value",
    }
    # autospec of configparser
    with mock.patch("configparser.ConfigParser", autospec=True) as config:
        _pc_wrapper.read_and_modify_config(**parameters)

        config.optionxform.assert_called_once
        config().read.assert_called_once_with(parameters["config_file_path"])
        config().has_section.assert_called_with(parameters["section"])
        config().set.assert_called_once_with(
            parameters["section"],
            parameters["parameter"],
            parameters["new_value"],
        )
        config().write.assert_called_once
