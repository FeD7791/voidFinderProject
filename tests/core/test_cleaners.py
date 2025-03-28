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

"""Test for cleaner methods."""

# =============================================================================
# IMPORTS
# =============================================================================

import os
import pathlib
import tempfile
from unittest import mock


import numpy as np

import pandas as pd

import pytest


from voidfindertk.core import cleaners
from voidfindertk.datasets import spherical_cloud


# =============================================================================
# TESTS
# =============================================================================


def _check_cbl_module():
    utils_path = pathlib.Path(os.path.dirname(os.path.abspath(__file__)))
    core_path = utils_path.parent.parent / "voidfindertk" / "core"
    return (core_path / "libcleaner.so").exists()


# @pytest.mark.skipif(
#     not _check_cbl_module(), reason="CBL module not available!"
# )
@pytest.mark.skip(reason="Requires CBL Installed")
def test_cbl_cleaner_interface(build_box_with_eq_voids):

    # Set initial variables
    rad = 30.0
    cloud = spherical_cloud.build_cloud()
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, rad=rad, delta=-0.9
    )

    # Build temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:

        # Get the path location of this Directory
        workdir = pathlib.Path(os.path.abspath(temp_dir))
        tracers_file_path = workdir / "tracers.txt"
        voids_file_path = workdir / "voids.txt"

        # Perfom file input creation.
        cleaners._save_xyz_tracers(box=box, path=tracers_file_path)
        cleaners._save_r_eff_center(
            centers=centers,
            r_eff=rad * np.ones(len(centers)),
            path=voids_file_path,
        )

        # Setting parameters to run the cleaner
        parameters = {
            "file_voids": voids_file_path,
            "file_tracers": tracers_file_path,
            "ratio": 1.5,
            "initial_radius": True,
            "delta_r_min": 10.0,
            "delta_r_max": 50.0,
            "threshold": 0.1,
            "output_path": workdir / "cleaned_catalogue.txt",
            "ol_crit": False,
            "rescale": True,
            "checkoverlap": True,
        }
        # Perform the cleaning using the CBL cleaner.
        cleaners._cbl_cleaner_interface(**parameters)


# @pytest.mark.skipif(
#     not _check_cbl_module(), reason="CBL module not available!"
# )
@pytest.mark.skip(reason="Requires CBL Installed")
def test_cbl_cleaner(build_box_with_eq_voids):
    # Setting parameters.
    rad = 30.0
    cloud = spherical_cloud.build_cloud()
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, rad=rad
    )

    # Build temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Get the path location of this Directory
        workdir = pathlib.Path(os.path.abspath(temp_dir))

        parameters1 = {
            "center": centers,
            "radius": rad * np.ones(len(centers)),
            "box": box,
            # temp directory
            "temporal_dir_path": workdir,
            "clean_directory": True,
            # cbl parameters
            "ratio": 1.5,
            "initial_radius": True,
            "delta_r_min": 10.0,
            "delta_r_max": 100.0,
            "threshold": 0.2,
            "ol_crit": "density_contrast",
            "rescale": True,
            "checkoverlap": True,
        }
        parameters2 = parameters1
        parameters2["ol_crit"] = "central_density"
        cbl_centers_1, cbl_radius_1 = cleaners._cbl_cleaner(**parameters1)
        cbl_centers_2, cbl_radius_2 = cleaners._cbl_cleaner(**parameters2)

    # Check we are getting reasonable results

    # Radius
    assert np.allclose(cbl_radius_1, rad, rtol=1.0)
    assert np.allclose(cbl_radius_2, rad, rtol=1.0)

    # Centers : cbl cleaners dismisses some of the voids as reals void, so it
    # ends cutting the catalogue a little bit


def test_overlap_cleaner():

    radius1 = np.array([30.0, 20.0, 40.0, 30.0, 20.0, 40.0, 5, 10, 2])
    centers1 = np.array(
        [
            [50.0, 50.0, 50.0],
            [60.0, 60.0, 60.0],
            [55.0, 55.0, 55.0],
            [150.0, 150.0, 150.0],
            [160.0, 160.0, 160.0],
            [155.0, 155.0, 155.0],
            [1500.0, 1500.0, 1500.0],
            [1501.0, 1501.0, 1501.0],
            [1502.0, 1502.0, 1502.0],
        ]
    )
    centers_new, radius_new = cleaners._overlap_cleaner(
        center=centers1, radius=radius1
    )

    assert all(
        np.all(
            centers_new
            == np.array(
                [
                    [55.0, 55.0, 55.0],
                    [155.0, 155.0, 155.0],
                    [1501.0, 1501.0, 1501.0],
                ]
            ),
            axis=1,
        )
    )
    assert np.all(radius_new == np.array([40.0, 40.0, 10.0]))

    # No overlapping case
    radius2 = np.array([30.0, 30.0, 30.0])
    centers2 = np.array(
        [
            [500.0, 500.0, 500.0],
            [600.0, 600.0, 600.0],
            [700.0, 700.0, 700.0],
        ]
    )
    centers_new2, radius_new2 = cleaners._overlap_cleaner(
        center=centers2, radius=radius2
    )
    assert all(np.all(centers_new2 == centers2, axis=1))
    assert np.all(radius_new2 == radius2)


def test_get_cleaner():

    with mock.patch(
        "voidfindertk.core.cleaners._overlap_cleaner"
    ) as mock_ov_cl:
        with mock.patch(
            "voidfindertk.core.cleaners._cbl_cleaner"
        ) as mock_cbl_cl:
            out1 = cleaners.get_cleaner(cleaner_method="overlap")
            out2 = cleaners.get_cleaner(cleaner_method="cbl")

    cleaner_method = "<fake_method>"

    with pytest.raises(
        ValueError, match=f"{cleaner_method} in not a valid cleaner method"
    ):
        cleaners.get_cleaner(cleaner_method=cleaner_method)

    assert mock_cbl_cl == out2
    assert mock_ov_cl == out1


def test_cbl_cleaner_mock_version():
    parameters = {
        "center": [1, 5, 9],
        "radius": [5, 6, 8],
        "box": "BOX",
        # temp directory
        "temporal_dir_path": ".",
        "clean_directory": True,
        # cbl parameters
        "ratio": 1.5,
        "initial_radius": True,
        "delta_r_min": 10.0,
        "delta_r_max": 100.0,
        "threshold": 0.2,
        "ol_crit": "density_contrast",
        "rescale": True,
        "checkoverlap": True,
    }
    directory_path = pathlib.Path("tempdir_path")
    cleaned_catalogue_path = (
        pathlib.Path(directory_path) / "cleaned_catalogue.txt"
    )

    input_centers_path = directory_path / "input_xyz_rad.txt"
    input_tracers_path = directory_path / "input_tracers.txt"
    with mock.patch("tempfile.mkdtemp", return_value=directory_path):
        with mock.patch("shutil.rmtree"):
            with mock.patch("numpy.linalg.norm"):
                with mock.patch(
                    "voidfindertk.core.cleaners._read_cleaned_catalogue",
                    new=mock.Mock(
                        return_value=[mock.MagicMock(), mock.MagicMock()]
                    ),
                ) as rcc_mock:
                    with mock.patch.multiple(
                        "voidfindertk.core.cleaners",
                        _save_xyz_tracers=mock.DEFAULT,
                        _save_r_eff_center=mock.DEFAULT,
                        _cbl_cleaner_interface=mock.DEFAULT,
                    ) as mocks:
                        cleaners._cbl_cleaner(**parameters)

    mocks["_save_xyz_tracers"].assert_called_once_with(
        box=parameters["box"], path=input_tracers_path
    )
    mocks["_save_r_eff_center"].assert_called_once_with(
        centers=parameters["center"],
        r_eff=parameters["radius"],
        path=input_centers_path,
    )

    mocks["_cbl_cleaner_interface"].assert_called_once_with(
        file_voids=input_centers_path,
        file_tracers=input_tracers_path,
        ratio=parameters["ratio"],
        initial_radius=parameters["initial_radius"],
        delta_r_min=parameters["delta_r_min"],
        delta_r_max=parameters["delta_r_max"],
        threshold=parameters["threshold"],
        output_path=cleaned_catalogue_path,
        # If ol_crit = "density_contrast" => this parameters is equal to
        # True inside function.
        ol_crit=True,
        rescale=parameters["rescale"],
        checkoverlap=parameters["checkoverlap"],
    )

    rcc_mock.assert_called_once_with(
        cleaned_catalogue_path=cleaned_catalogue_path
    )


def test_cbl_cleaner_interface_mock_version():
    # Mock the ctypes.CDLL to return a mock object
    with mock.patch("ctypes.CDLL") as mock_cdll:
        # Create a mock for the 'clibrary' object returned by ctypes.CDLL
        mock_clibrary = mock.Mock()

        # Assign the mock object to be returned when ctypes.CDLL is called
        mock_cdll.return_value = mock_clibrary

        # Mock the 'process_catalogues' function inside the clibrary mock
        mock_clibrary.process_catalogues = mock.Mock()

        # Test data for the function call
        parameters = {
            "file_voids": "voids_file.txt",
            "file_tracers": "tracers_file.txt",
            "ratio": 0.5,
            "initial_radius": True,
            "delta_r_min": 10.0,
            "delta_r_max": 100.0,
            "threshold": 0.2,
            "output_path": "output_catalogue.txt",
            "ol_crit": True,
            "rescale": True,
            "checkoverlap": True,
        }

        # Call the function being tested
        cleaners._cbl_cleaner_interface(**parameters)

        # Assert that the 'process_catalogues' method was called with the
        # correct arguments
        mock_clibrary.process_catalogues.assert_called_once_with(
            parameters["file_voids"].encode("utf-8"),
            parameters["file_tracers"].encode("utf-8"),
            parameters["ratio"],
            parameters["initial_radius"],
            mock.ANY,
            mock.ANY,
            parameters["threshold"],
            parameters["output_path"].encode("utf-8"),
            parameters["ol_crit"],
            parameters["rescale"],
            parameters["checkoverlap"],
        )


def test_save_xyz_tracers():
    mock_box = mock.MagicMock()
    with mock.patch("numpy.column_stack"):
        with mock.patch("pandas.DataFrame") as mock_df:
            cleaners._save_xyz_tracers(box=mock_box, path=pathlib.Path("path"))
    mock_df.assert_called_once()


def test_save_r_eff_center():
    centers = np.ones((15, 3))
    r_eff = 30.0 * np.ones((15,))
    with tempfile.TemporaryDirectory() as tempdir:
        path = pathlib.Path(tempdir) / "file.txt"
        cleaners._save_r_eff_center(centers=centers, path=path, r_eff=r_eff)
        # Now we will retrieve this variables
        centers_r_eff = pd.read_csv(
            path, delim_whitespace=True, names=["x", "y", "z", "r"]
        )

    assert centers_r_eff.shape == (15, 4)
    assert np.all(np.array(centers_r_eff[["x", "y", "z"]]) == centers)
    assert np.all(np.array(centers_r_eff["r"]) == r_eff)
