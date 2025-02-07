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
    parameters = {
        "mpi_flags": "mpi_flags",
        "bin_path": "bin_path",
        "conf_file_path": "conf_file_path",
        "work_dir_path": ".",
    }

    popcorn_mock = mock.MagicMock()
    with mock.patch("sh.Command", return_value=popcorn_mock) as shcommand:
        _pc_wrapper.popcorn_void_finder(**parameters)

    shcommand.assert_called_once_with(
        "popcorn", search_paths=[parameters["bin_path"]]
    )
    popcorn_mock.assert_called_once_with(
        "config=" + str(parameters["conf_file_path"])
    )


def test_compute_intersects():
    parameters = {
        "bin_path": "bin_path",
        "conf_file_path": "conf_file_path",
        "work_dir_path": ".",
    }
    compute_intersects_mock = mock.MagicMock()

    with mock.patch(
        "sh.Command", return_value=compute_intersects_mock
    ) as shcommand:
        _pc_wrapper.compute_intersects(**parameters)

    shcommand.assert_called_once_with(
        "compute_intersecs", search_paths=[parameters["bin_path"]]
    )

    compute_intersects_mock.assert_called_once_with(
        "config=" + parameters["conf_file_path"]
    )


def test_clean_duplicates():
    parameters = {
        "bin_path": "bin_path",
        "conf_file_path": "conf_file_path",
        "work_dir_path": ".",
    }
    clean_duplicates_mock = mock.MagicMock()

    with mock.patch(
        "sh.Command", return_value=clean_duplicates_mock
    ) as shcommand:
        _pc_wrapper.clean_duplicates(**parameters)

    shcommand.assert_called_once_with(
        "clean_duplicates", search_paths=[parameters["bin_path"]]
    )

    clean_duplicates_mock.assert_called_once_with(
        "config=" + parameters["conf_file_path"]
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
