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

import numpy as np

from voidfindertk.core import cleaners
from voidfindertk.datasets import spherical_cloud

# =============================================================================
# TESTS
# =============================================================================


def test_cbl_cleaner_interface(build_box_with_eq_voids):

    # Set initial variables
    rad = 30.0
    cloud = spherical_cloud.build_cloud()
    box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
        cloud=cloud, rad=rad
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
            "threshold": 0.3,
            "output_path": workdir / "cleaned_catalogue.txt",
            "ol_crit": True,
            "rescale": True,
            "checkoverlap": True,
        }
        # Perform the cleaning using the CBL cleaner.
        cleaners._cbl_cleaner_interface(**parameters)


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
            "ratio": 1.0,
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
