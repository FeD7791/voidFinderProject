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

"""Test PopCorn."""

# =============================================================================
# IMPORTS
# =============================================================================

import pathlib
from unittest import mock


import pytest


from voidfindertk.datasets import spherical_cloud
from voidfindertk.popcorn import PopCorn
from voidfindertk.svf_popcorn import SVFPopCorn

# =============================================================================
# TEST
# =============================================================================


def test_popcorn_model_find(mkbox):

    params = {
        "auxfiles": "true",
        "boxsize": 1000.0,
        "densth": -0.9,
        "minradius": 5,
        "maxradius": 100,
        "svf_path": pathlib.Path("."),
        "workdir": pathlib.Path("."),
        "workdir_clean": False,
        "shot_noise_threshold": 20,
    }
    box_ = mkbox(seed=42)

    with mock.patch.object(SVFPopCorn, "model_find") as svf_popcorn_mock:
        model_find_parameters = PopCorn(**params).model_find(box=box_)

    svf_popcorn_mock.assert_called_once()
    svf_popcorn_mock.assert_called_once_with(box_)
    model_find_parameters.__setitem__.assert_called_once_with(
        "build_popcorn", True
    )


def test_popcorn_build_voids_true(svf_popcorn_paths_and_names, mkbox):
    params = {
        "auxfiles": "true",
        "boxsize": 1000.0,
        "densth": -0.9,
        "minradius": 5,
        "maxradius": 100,
        "svf_path": pathlib.Path("."),
        "workdir": pathlib.Path("."),
        "workdir_clean": False,
        "shot_noise_threshold": 20,
    }
    box_ = mkbox(seed=42)
    model_find_parameters = {
        "build_popcorn": True,
        "box": box_,
        "run_work_dir": ".",
    }
    run_work_dir = pathlib.Path(".")
    pn = svf_popcorn_paths_and_names()

    # Useful Mocks
    voids = mock.MagicMock()
    spheres = mock.MagicMock()

    with mock.patch.object(SVFPopCorn, "model_find"):
        with mock.patch.object(SVFPopCorn, "build_voids") as svf_build_voids:
            svf_build_voids.return_value = [
                mock.MagicMock(),
                mock.MagicMock(),
                {"files_directory_path": "."},
            ]
            with mock.patch.multiple(
                "voidfindertk.popcorn._pc_wrapper",
                read_and_modify_config=mock.DEFAULT,
                popcorn_void_finder=mock.DEFAULT,
                compute_intersects=mock.DEFAULT,
                clean_duplicates=mock.DEFAULT,
            ) as _pc_wrapper_mock:

                with mock.patch(
                    "voidfindertk.popcorn._pc_postprocessing.get_properties",
                    return_value=[voids, spheres],
                ) as get_properties_mock:

                    # Run the function
                    tracers, popcorn_centers, extra = PopCorn(
                        **params
                    ).build_voids(model_find_parameters=model_find_parameters)

    # Testing assertions for build_popcorn = True
    _pc_wrapper_mock["read_and_modify_config"].assert_called_once_with(
        config_file_path=run_work_dir / pn.CONFIG,
        section="INPUT_PARAMS",
        parameter="MINRADIUS",
        new_value=str(params["shot_noise_threshold"]),
    )
    _pc_wrapper_mock["popcorn_void_finder"].assert_called_once_with(
        bin_path=pn.SVF,
        conf_file_path=run_work_dir / pn.CONFIG,
        work_dir_path=run_work_dir,
    )
    _pc_wrapper_mock["compute_intersects"].assert_called_once_with(
        bin_path=pn.SVF,
        conf_file_path=run_work_dir / pn.CONFIG,
        work_dir_path=run_work_dir,
    )
    _pc_wrapper_mock["clean_duplicates"].assert_called_once_with(
        bin_path=pn.SVF,
        conf_file_path=run_work_dir / pn.CONFIG,
        work_dir_path=run_work_dir,
    )
    get_properties_mock.assert_called_once_with(
        filename=run_work_dir / pn.POPFILE
    )

    voids.sort_values.assert_called_once_with(by=["id"], inplace=True)
    spheres.sort_values.assert_called_once_with(by=["id"], inplace=True)
    assert voids == extra["voids"]
    assert spheres == extra["spheres"]
    assert tracers == voids.tracers
    assert run_work_dir == extra["files_directory_path"]


def test_popcorn_build_voids_false(mkbox):
    params = {
        "auxfiles": "true",
        "boxsize": 1000.0,
        "densth": -0.9,
        "minradius": 5,
        "maxradius": 100,
        "svf_path": pathlib.Path("."),
        "workdir": pathlib.Path("."),
        "workdir_clean": False,
        "shot_noise_threshold": 20,
    }
    box_ = mkbox(seed=42)
    model_find_parameters = {
        "build_popcorn": False,
        "box": box_,
        "run_work_dir": ".",
    }
    run_work_dir = pathlib.Path(model_find_parameters["run_work_dir"])

    # Useful Mocks
    voids = mock.MagicMock()
    spheres = mock.MagicMock()

    # Testing for build_popcorn: False
    model_find_parameters = {
        "build_popcorn": False,
        "box": box_,
        "run_work_dir": ".",
    }
    with mock.patch(
        "voidfindertk.popcorn._pc_postprocessing.get_properties",
        return_value=[voids, spheres],
    ):

        tracers, popcorn_centers, extra = PopCorn(**params).build_voids(
            model_find_parameters=model_find_parameters
        )

    voids.sort_values.assert_called_once_with(by=["id"], inplace=True)
    spheres.sort_values.assert_called_once_with(by=["id"], inplace=True)
    assert voids == extra["voids"]
    assert spheres == extra["spheres"]
    assert tracers == voids.tracers
    assert run_work_dir == extra["files_directory_path"]


# @pytest.mark.skipif(
#     not (pathlib.Path(SETTINGS.popcorn_path) / "svf").exists(),
#     reason="POPCORN not available!",
# )
@pytest.mark.skip(reason="Requires POPCORN Installed")
def test_popcorn_working_example(build_box_with_eq_voids):

    delta = -0.9
    cloud = spherical_cloud.build_cloud(lmin=0, lmax=1000)
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, delta=delta, rad=30
    )
    model = PopCorn(
        workdir_clean=True,
        boxsize=1000,
        densth=delta,
        minradius=5,
        maxradius=100,
        shot_noise_threshold=15,
        cores=4,
    )

    model.find(box=box)


@pytest.mark.skip(reason="Requires POPCORN Installed")
def test_popcorn_box_tracers_mapping(build_box_with_eq_voids):

    delta = -0.8
    cloud = spherical_cloud.build_cloud(lmin=0, lmax=1000, n_points=100**3)
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, delta=delta, rad=30, log_mass_min=10, log_mass_max=50
    )

    model = PopCorn(
        workdir_clean=True,
        boxsize=1000,
        densth=delta,
        minradius=5,
        maxradius=100,
        shot_noise_threshold=15,
    )

    void = model.find(box=box)
    df = void.extra_.voids
    df["tracers2"] = void.tracers
    assert all(df.tracers == df.tracers2)
