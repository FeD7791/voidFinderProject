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

import ctypes
import pathlib
from unittest import mock


from voidfindertk.zobov import _zb_wrapper

# =============================================================================
# TESTS
# =============================================================================


def test_run_vozinit(load_mock_data):

    params = {
        "input_file_path":pathlib.Path("."),
        "buffer_size":2,
        "box_size":1000,
        "number_of_divisions":2,
        "executable_name":"vozinit",
        "work_dir_path":pathlib.Path("./workdir"),
        "vozinit_dir_path": pathlib.Path("./vozinit")
    }


    # mock the class and run
    with mock.patch("sh.Command") as shcommand:
        _zb_wrapper.run_vozinit(**params)

    shcommand(
        "vozinit", search_paths=[params["vozinit_dir_path"]]
    ).assert_called_once_with(
        str(params["input_file_path"]),
        str(params["buffer_size"]),
        str(params["box_size"]),
        str(params["number_of_divisions"]),
        str(params["executable_name"]),
        _cwd=params["work_dir_path"],
        _out=shcommand("vozinit").call_args[1]["_out"],
        _err_to_out=True,
    )


def test_run_voz_step():

    params = {
        "preprocess_dir_path": pathlib.Path("./path/dir"),
        "executable_name": "<EXE_NAME>",
        "work_dir_path": pathlib.Path("."),
        "voz_executables_path": pathlib.Path("./exe/path"),
    }
    # Expected calls
    expected_call_1 = mock.call(
        params["voz_executables_path"] / "voz1b1", params["work_dir_path"]
    )
    expected_call_2 = mock.call(
        params["voz_executables_path"] / "voztie", params["work_dir_path"]
    )

    mock_preprocess = mock.MagicMock()

    with mock.patch("voidfindertk.zobov._zb_wrapper._move_inputs") as mock_mi:
        mock_mi.return_value = mock_preprocess

        with mock.patch("sh.Command") as mock_command:
            mock_command.return_value = mock_preprocess
            _zb_wrapper.run_voz_step(**params)

    mock_command.assert_called_once_with(
        str(params["preprocess_dir_path"] / f"scr{params['executable_name']}")
    )
    mock_preprocess.assert_called_once()
    mock_mi.assert_has_calls(
        [expected_call_1, expected_call_2], any_order=False
    )
    assert mock_mi.call_count == 2


def test_run_voz1b1():
    """Tests run_voz1b1"""
    params = {
        "input_file_path": "input_file_path",
        "buffer_size": "buffer_size",
        "box_size": "box_size",
        "executable_name": "executable_name",
        "number_of_divisions": "number_of_divisions",
        "binary_division": [0, 0, 0],
        "voz1b1_dir_path": pathlib.Path("voz1b1_dir_path"),
        "work_dir_path": ".",
    }
    args = [
        str(params["input_file_path"]),
        str(params["buffer_size"]),
        str(params["box_size"]),
        str(params["executable_name"]),
        str(params["number_of_divisions"]),
        str(params["binary_division"][0]),
        str(params["binary_division"][1]),
        str(params["binary_division"][2]),
    ]

    with mock.patch("sh.Command") as shcommand:

        _zb_wrapper.run_voz1b1(**params)

    shcommand(
        "voz1b1", search_paths=[params["voz1b1_dir_path"]]
    ).assert_called_once_with(
        *args,
        _cwd=params["work_dir_path"],
        _out=shcommand("mpirun").call_args[1]["_out"],
        _err_to_out=True,
    )


def test_run_voztie():
    """Test run_voztie"""
    params = {
        "number_of_divisions": "number_of_divisions",
        "executable_name": "executable_name",
        "voztie_dir_path": pathlib.Path("voztie_dir_path"),
        "work_dir_path": ".",
    }

    with mock.patch("sh.Command") as shcommand:
        _zb_wrapper.run_voztie(**params)

    shcommand(
        "voztie", search_paths=[params["voztie_dir_path"]]
    ).assert_called_once_with(
        params["number_of_divisions"],
        params["executable_name"],
        _cwd=params["work_dir_path"],
        _out=shcommand("mpirun").call_args[1]["_out"],
        _err_to_out=True,
    )


def test_run_jozov():
    params = {
        "jozov_dir_path": pathlib.Path("jozov_dir_path"),
        "executable_name": "executable_name",
        "output_name_particles_in_zones": "output_name_particles_in_zones",
        "output_name_zones_in_void": "output_name_zones_in_void",
        "output_name_text_file": "output_name_text_file",
        "density_threshold": "density_threshold",
        "work_dir_path": ".",
    }
    args = (
        f"adj{params['executable_name']}.dat",
        f"vol{params['executable_name']}.dat",
        f"{params['output_name_particles_in_zones']}.dat",
        f"{params['output_name_zones_in_void']}.dat",
        f"{params['output_name_text_file']}.dat",
        f"{params['density_threshold']}",
    )

    with mock.patch("sh.Command") as shcommand:
        _zb_wrapper.run_jozov(**params)

    shcommand(
        "jozov", search_paths=[params["jozov_dir_path"]]
    ).assert_called_once_with(
        *args,
        _cwd=params["work_dir_path"],
        _out=shcommand("jozov").call_args[1]["_out"],
        _err_to_out=True,
    )


def test_write_input():
    mock_box = mock.MagicMock()
    params = {
        "box": mock_box,
        "path_executable": "path_executable",
        "raw_file_path": "raw_file_path",
        "txt_file_path": "txt_file_path",
    }
    clibrary = mock.MagicMock()

    with mock.patch("ctypes.CDLL") as mock_cdll:
        mock_cdll.return_value = clibrary
        _zb_wrapper.write_input(**params)

    mock_cdll.assert_called_once_with(
        str(params["path_executable"]), mode=ctypes.RTLD_GLOBAL
    )
    clibrary.c_binary_writter.assert_called_once_with(
        mock_box.arr_.x,
        mock_box.arr_.y,
        mock_box.arr_.z,
        len(mock_box),
        str(params["raw_file_path"]).encode("utf-8"),
        str(params["txt_file_path"]).encode("utf-8"),
    )
