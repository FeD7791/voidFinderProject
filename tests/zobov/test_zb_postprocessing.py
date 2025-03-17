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

"""
Test ZOBOV posprocessing modules

"""

# =============================================================================
# IMPORTS
# =============================================================================

import ctypes
import os
import pathlib
import tempfile
from unittest import mock


import numpy as np

import pandas as pd

import pytest


from voidfindertk.zobov import _zb_postprocessing


# =============================================================================
# TESTS
# =============================================================================


def test_parse_zones_in_void_output():
    parameters = {
        "executable_path": "executable_path",
        "input_file_path": "input_file_path",
        "output_file_path": "output_file_path",
    }
    clibrary = mock.MagicMock()

    with mock.patch("ctypes.CDLL") as mock_cdll:
        mock_cdll.return_value = clibrary
        _zb_postprocessing.parse_zones_in_void_output(**parameters)
    mock_cdll.assert_called_once_with(
        str(parameters["executable_path"]), mode=ctypes.RTLD_GLOBAL
    )
    clibrary.process_files.assert_called_once_with(
        str(parameters["input_file_path"]).encode(),
        str(parameters["output_file_path"]).encode(),
    )


def test_parse_tracers_in_zones_output():
    parameters = {
        "executable_path": "executable_path",
        "input_file_path": "input_file_path",
        "output_file_path": "output_file_path",
    }
    clibrary = mock.MagicMock()

    with mock.patch("ctypes.CDLL") as mock_cdll:
        mock_cdll.return_value = clibrary
        _zb_postprocessing.parse_tracers_in_zones_output(**parameters)
    mock_cdll.assert_called_once_with(
        str(parameters["executable_path"]), mode=ctypes.RTLD_GLOBAL
    )
    clibrary.get_tracers_in_zones.assert_called_once_with(
        str(parameters["input_file_path"]).encode(),
        str(parameters["output_file_path"]).encode(),
    )


def test_get_tracers_in_zones():

    tests_path = pathlib.Path(os.path.abspath(__file__)).parent.parent
    tracers_vs_zones_file = (
        tests_path / "mock_data" / "zobov" / "part_vs_zone_ascii.txt"
    )
    with tempfile.TemporaryDirectory():

        out = pd.read_csv(
            tests_path / "mock_data" / "zobov" / "output_txt.dat",
            header=1,
            delim_whitespace=True,
        )
        # Tracers that are in each zone
        tracers = _zb_postprocessing._get_tracers_in_zones(
            tracers_in_zones_path=tracers_vs_zones_file
        )

    # This have the same order:
    core_particle = np.array([t[0] for t in tracers])
    out_core_particle = np.array(out["CoreParticle"])

    indx = np.argsort(out_core_particle)

    # the number of tracers in tracers should be the same as out['Zone#Part']
    # when they are in the right order
    number_of_tracers = np.array(list(map(len, tracers)))

    assert np.all(core_particle == out_core_particle[indx])
    assert np.all(number_of_tracers == np.array(out["Zone#Part"][indx]))


def test_get_zones_in_void():

    tests_path = pathlib.Path(os.path.abspath(__file__)).parent.parent
    zones_vs_voids_file = (
        tests_path / "mock_data" / "zobov" / "zones_vs_voids_ascii.txt"
    )
    with tempfile.TemporaryDirectory():

        out = pd.read_csv(
            tests_path / "mock_data" / "zobov" / "output_txt.dat",
            header=1,
            delim_whitespace=True,
        )

        void_number = np.array(out["Void#"])
        zones = _zb_postprocessing._get_zones_in_void(
            zones_in_void_file_path=zones_vs_voids_file
        )
    zone_number = np.array([z[0] for z in zones]) + 1
    number_of_zones = np.array(out["Void#Zones"])
    # You can think of 'two' columns 1: #Zone 2: then #Zone + other zones
    number_of_zones_from_module = np.array(list(map(len, zones))) - 1

    assert np.all(void_number == zone_number)
    assert np.all(
        np.sort(number_of_zones) == np.sort(number_of_zones_from_module)
    )


def test_get_tracers_in_voids():

    tests_path = pathlib.Path(os.path.abspath(__file__)).parent.parent
    zones_vs_voids_file = (
        tests_path / "mock_data" / "zobov" / "zones_vs_voids_ascii.txt"
    )
    tracers_vs_zones_file = (
        tests_path / "mock_data" / "zobov" / "part_vs_zone_ascii.txt"
    )

    with tempfile.TemporaryDirectory():

        out = pd.read_csv(
            tests_path / "mock_data" / "zobov" / "output_txt.dat",
            header=1,
            delim_whitespace=True,
        )
        (properties_df, tracers_in_void) = (
            _zb_postprocessing.get_tracers_in_voids(
                properties_dataframe=out,
                tracers_in_zones_path=tracers_vs_zones_file,
                zones_in_void_path=zones_vs_voids_file,
            )
        )
    number_tracers_in_void = np.array(list(map(len, tracers_in_void)))
    df_n_tracers_in_void = np.array(properties_df["Void#Part"])

    assert np.all(number_tracers_in_void == df_n_tracers_in_void)


def test_get_center_method():

    with mock.patch(
        "voidfindertk.zobov._zb_postprocessing._centers_barycentre_method"
    ) as mock_centers_barycentre_method:
        with mock.patch(
            "".join(
                [
                    "voidfindertk.zobov._zb_postprocessing.",
                    "_center_core_particle_method",
                ]
            )
        ) as mock_center_core_particle_method:
            output_1 = _zb_postprocessing.get_center_method(
                center_method="barycentre"
            )
            output_2 = _zb_postprocessing.get_center_method(
                center_method="core_particle"
            )
    with pytest.raises(
        ValueError, match="This center_method is not available!"
    ):
        _zb_postprocessing.get_center_method(center_method="undefined_method")
    assert mock_centers_barycentre_method == output_1
    assert mock_center_core_particle_method == output_2


def test_center_core_particle_method(mkbox):
    box = mkbox()
    properties_df = {"CoreParticle": pd.Series([0, 1, 5])}

    xyz = _zb_postprocessing._center_core_particle_method(
        properties_df=properties_df, box=box
    )

    assert xyz[0][0] == box.arr_.x[0]
    assert xyz[1][0] == box.arr_.x[1]
    assert xyz[2][0] == box.arr_.x[5]

    assert xyz[0][1] == box.arr_.y[0]
    assert xyz[1][1] == box.arr_.y[1]
    assert xyz[2][1] == box.arr_.y[5]

    assert xyz[0][2] == box.arr_.z[0]
    assert xyz[1][2] == box.arr_.z[1]
    assert xyz[2][2] == box.arr_.z[5]


def test_read_volume_file():
    # The volumes are from the voronoi cells, there were 700070 tracers
    # then volumes should have 700070 elements.
    path = pathlib.Path(os.path.abspath(__file__)).parent.parent

    volumes = _zb_postprocessing._read_volume_file(
        filename=path / "mock_data" / "zobov" / "voloutput_vozinit.dat"
    )
    assert volumes.shape == (700070,)


def test_get_tracers_xyz(mkbox):
    box = mkbox()
    assert _zb_postprocessing._get_tracers_xyz(box=box).shape == (len(box), 3)


# def test_centers_barycentre_method():

#     with tempfile.TemporaryDirectory() as tempdir:
#         model = zobov.ZobovVF(
#             workdir=workdir_path,
#             box_size=1000,
# )
