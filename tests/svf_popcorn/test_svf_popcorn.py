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

"""Test for voidfindertk.svf_popcorn.SVFPopCorn."""

# =============================================================================
# IMPORTS
# =============================================================================

import pathlib
from unittest import mock


import pytest


from voidfindertk.datasets import spherical_cloud
from voidfindertk.svf_popcorn import SVFPopCorn


# =============================================================================
# TESTS
# =============================================================================


@pytest.mark.skip(reason="Requires POPCORN Installed")
@pytest.mark.parametrize("cores", [8])
def test_svf_popcorn_working_example(cores, build_box_with_eq_voids):
    """
    Tests svfpopcorn works with a real run
    """
    delta = -0.9
    cloud = spherical_cloud.build_cloud(lmin=0, lmax=1000)
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, delta=delta, rad=30
    )
    model = SVFPopCorn(
        workdir_clean=True,
        boxsize=1000,
        densth=delta,
        minradius=5,
        maxradius=100,
        cores=cores,
    )

    model.find(box=box)


def test_svfpopcorn(svf_popcorn_paths_and_names, mkbox):

    params = {
        "auxfiles": "true",
        "boxsize": 1000.0,
        "densth": -0.9,
        "minradius": 5,
        "maxradius": 100,
        "svf_path": pathlib.Path("."),
        "workdir": pathlib.Path("."),
        "workdir_clean": False,
        "cores": 2,
    }
    box_ = mkbox(seed=42)
    pn = svf_popcorn_paths_and_names()

    with mock.patch(
        "voidfindertk.utils.make_workdir.create_run_work_dir",
        return_value=mock.MagicMock(),
    ) as tempdir:
        with mock.patch.multiple(
            "voidfindertk.svf_popcorn._svf_pc_wrapper",
            config_file_maker=mock.DEFAULT,
            popcorn_svf_input_data_builder=mock.DEFAULT,
            spherical_popcorn_void_finder=mock.DEFAULT,
        ) as svf_pc_wrapper_mocks:
            SVFPopCorn(**params).model_find(box=box_)

            with mock.patch.multiple(
                "voidfindertk.svf_popcorn._svf_pc_postprocessing",
                get_void_properties=mock.DEFAULT,
                get_tracers_in_voids=mock.DEFAULT,
            ) as svf_pc_postprocessing_mocks:
                tinv, xyz_properties, extra = SVFPopCorn(**params).build_voids(
                    model_find_parameters={
                        "run_work_dir": pathlib.Path("./run_work_dir"),
                        "box": box_,
                    }
                )

    svf_pc_wrapper_mocks["config_file_maker"].assert_called_once_with(
        trsfile=str(tempdir(workdir_path=params["workdir"]) / pn.TRSFILE),
        filefmt="ASCII",
        num_file=str(1),
        sphfile=str(tempdir(workdir_path=params["workdir"]) / pn.SPHFILE),
        popfile=str(tempdir(workdir_path=params["workdir"]) / pn.POPFILE),
        auxfiles=str(params["auxfiles"]),
        rawpopfile=str(
            tempdir(workdir_path=params["workdir"]) / pn.RAWPOPFILE
        ),
        pairsfile=str(tempdir(workdir_path=params["workdir"]) / pn.PAIRSFILE),
        boxsize=str(params["boxsize"]),
        densth=str(params["densth"]),
        minradius=str(params["minradius"]),
        maxradius=str(params["maxradius"]),
        massmin=str(0),  # This is always the value for this parameter.
        eps=str(1e-5),
        path=str(tempdir(workdir_path=params["workdir"]) / pn.CONFIG),
    )

    svf_pc_wrapper_mocks[
        "popcorn_svf_input_data_builder"
    ].assert_called_once_with(
        box=box_,
        file_path=str(tempdir(workdir_path=params["workdir"]) / pn.TRSFILE),
    )
    svf_pc_wrapper_mocks[
        "spherical_popcorn_void_finder"
    ].assert_called_once_with(
        bin_path=pn.SVF,
        conf_file_path=tempdir(workdir_path=params["workdir"]) / pn.CONFIG,
        work_dir_path=tempdir(workdir_path=params["workdir"]),
        cores=params["cores"],
    )
    svf_pc_postprocessing_mocks["get_void_properties"].assert_called_once_with(
        popcorn_output_file_path=str(
            pathlib.Path("./run_work_dir") / pn.SPHFILE
        )
    )
    svf_pc_postprocessing_mocks[
        "get_tracers_in_voids"
    ].assert_called_once_with(
        box=box_,
        popcorn_output_file_path=str(
            pathlib.Path("./run_work_dir") / pn.SPHFILE
        ),
    )
