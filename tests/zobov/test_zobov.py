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

"""Tests for voidfindertk.zobov.ZobovVF."""

# =============================================================================
# IMPORTS
# =============================================================================

import pathlib
from unittest import mock


import numpy as np

import pytest


from voidfindertk.datasets import spherical_cloud
from voidfindertk.zobov import ZobovVF


def test_zobovvf(zobov_paths_and_names, mkbox):
    zpn = zobov_paths_and_names()
    params = {
        "buffer_size": 0.08,
        "box_size": 500,
        "number_of_divisions": 2,
        "density_threshold": 0.2,
        "zobov_path": pathlib.Path("."),
        "workdir": pathlib.Path("."),
        "workdir_clean": False,
        "dtype": np.float32,
    }
    run_workdir = pathlib.Path("run_workdir")

    box_ = mkbox(seed=42)
    mock_df = mock.MagicMock()
    mock_tinv = mock.MagicMock()
    properties = {"CoreParticle": [0, 1, 2]}
    extra_ = {
        "zobov_path": params["zobov_path"],
        "properties": properties,
        "files_directory_path": run_workdir,
    }

    with mock.patch("tempfile.mkdtemp") as mock_mkdetemp:
        mock_mkdetemp.return_value = pathlib.Path("run_workdir")
        with mock.patch("pandas.read_csv") as mock_pd_read_csv:
            mock_pd_read_csv.return_value = mock_df
            with mock.patch.multiple(
                "voidfindertk.zobov._zb_wrapper",
                # if don't use mock.DEFAULT dict is empty
                write_input=mock.DEFAULT,
                run_vozinit=mock.DEFAULT,
                run_voz_step=mock.DEFAULT,
                run_jozov=mock.DEFAULT,
            ) as mocks:
                ZobovVF(**params).model_find(box=box_)
                with mock.patch.multiple(
                    "voidfindertk.zobov._zb_postprocessing",
                    parse_tracers_in_zones_output=mock.DEFAULT,
                    parse_zones_in_void_output=mock.DEFAULT,
                    get_tracers_in_voids=mock.DEFAULT,
                    get_center_method=mock.DEFAULT,
                ) as mocks2:
                    mocks2["get_tracers_in_voids"].return_value = (
                        # properties, tinv
                        properties,
                        mock_tinv,
                    )
                    tinv, centers, extra = ZobovVF(**params).build_voids(
                        model_find_parameters={
                            "run_work_dir": run_workdir,
                            "box": box_,
                            "center_method": "barycentre",
                        }
                    )

    mocks["write_input"].assert_called_once_with(
        box=box_,
        path_executable=zpn.CURRENT_FILE_PATH / zpn.ZOBOV_LOADER_BIN,
        raw_file_path=run_workdir / zpn.TRACERS_RAW,
        txt_file_path=run_workdir / zpn.TRACERS_TXT,
    )
    mocks["run_vozinit"].assert_called_once_with(
        vozinit_dir_path=params["zobov_path"] / "src",
        input_file_path=run_workdir / zpn.TRACERS_RAW,
        buffer_size=params["buffer_size"],
        box_size=params["box_size"],
        number_of_divisions=params["number_of_divisions"],
        executable_name=zpn.OUTPUT_VOZINIT,
        work_dir_path=run_workdir,
    )
    mocks["run_voz_step"].assert_called_once_with(
        preprocess_dir_path=run_workdir,
        executable_name=zpn.OUTPUT_VOZINIT,
        work_dir_path=run_workdir,
        voz_executables_path=params["zobov_path"] / "src",
    )
    mocks["run_jozov"].assert_called_once_with(
        jozov_dir_path=params["zobov_path"] / "src",
        executable_name="output_vozinit",
        output_name_particles_in_zones=zpn.PARTICLES_IN_ZONES,
        output_name_zones_in_void=zpn.ZONES_IN_VOID,
        output_name_text_file=zpn.OUTPUT_JOZOV_VOIDS,
        density_threshold=params["density_threshold"],
        work_dir_path=run_workdir,
    )

    mocks2["parse_tracers_in_zones_output"].assert_called_once_with(
        executable_path=zpn.CURRENT_FILE_PATH / zpn.TRACERS_IN_ZONES_BIN,
        input_file_path=run_workdir / zpn.PARTICLES_VS_ZONES_RAW,
        output_file_path=run_workdir / zpn.PARTICLES_VS_ZONES_ASCII,
    )

    mocks2["parse_zones_in_void_output"].assert_called_once_with(
        executable_path=zpn.CURRENT_FILE_PATH / zpn.ZONES_IN_VOIDS_BIN,
        input_file_path=run_workdir / zpn.ZONES_VS_VOID_RAW,
        output_file_path=run_workdir / zpn.ZONES_VS_VOID_ASCII,
    )

    mocks2["get_tracers_in_voids"].assert_called_once_with(
        properties_dataframe=mock_df,
        tracers_in_zones_path=(run_workdir / zpn.PARTICLES_VS_ZONES_ASCII),
        zones_in_void_path=run_workdir / zpn.ZONES_VS_VOID_ASCII,
    )

    assert tinv == tuple(mock_tinv)
    assert np.all(
        mocks2["get_center_method"]("core_particle").return_value == centers
    )
    assert extra == extra_


@pytest.mark.skip(reason="Requires ZOBOV Installed")
def test_zobov_working_example(build_box_with_eq_voids):
    delta = -0.8
    cloud = spherical_cloud.build_cloud(lmin=0, lmax=1000)
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, delta=delta, rad=30
    )
    model = ZobovVF(
        workdir_clean=True,
        box_size=box.size(),
    )
    model.find(box=box)
