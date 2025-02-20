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

"""Module that holds functions and methods that are used to run the Popcorn \
void finder."""

# =============================================================================
# IMPORTS
# =============================================================================

import pathlib

import attr

import numpy as np

from . import _pc_postprocessing, _pc_wrapper
from ..svf_popcorn import FileNames, Paths, SVFPopCorn


# =============================================================================
# VFINDER
# =============================================================================


@attr.define
class PopCorn(SVFPopCorn):
    """
    A class to represent a popcorn void finder which processes spatial data.

    Attributes
    ----------
    _shot_noise_threshold : int
        The threshold value for shot noise, default is 20.

    Methods
    -------
    shot_noise_threshold
        Property that returns the shot noise threshold.
    build_voids(model_find_parameters)
        Builds voids in the model using specified parameters and updates
        configuration as necessary.
    """

    _shot_noise_threshold = attr.field(
        default=20, alias="shot_noise_threshold"
    )

    @property
    def shot_noise_threshold(self):
        """Return the shot noise threshold value."""
        return self._shot_noise_threshold

    def model_find(self, box):
        """Performs PopCorn Search."""
        model_find_parameters = super().model_find(box)
        model_find_parameters["build_popcorn"] = True
        return model_find_parameters

    def build_voids(self, model_find_parameters):
        """
        Builds voids in the model based on provided parameters.

        Parameters
        ----------
        model_find_parameters : dict
            A dictionary of parameters used to locate model features.

        Returns
        -------
        tuple
            A tuple containing tracers, centers, and additional information
            regarding the files directory and effective radii.
        """
        if model_find_parameters["build_popcorn"]:
            tracers_in_voids, centers, extra = super().build_voids(
                model_find_parameters
            )
            run_work_dir = pathlib.Path(extra["files_directory_path"])

            # Before continuing minradius must be re-configured
            _pc_wrapper.read_and_modify_config(
                config_file_path=run_work_dir / FileNames.CONFIG,
                section="INPUT_PARAMS",
                parameter="MINRADIUS",
                new_value=str(self._shot_noise_threshold),
            )
            _pc_wrapper.popcorn_void_finder(
                bin_path=Paths.SVF,
                conf_file_path=run_work_dir / FileNames.CONFIG,
                work_dir_path=run_work_dir,
            )
            _pc_wrapper.compute_intersects(
                bin_path=Paths.SVF,
                conf_file_path=run_work_dir / FileNames.CONFIG,
                work_dir_path=run_work_dir,
            )
            _pc_wrapper.clean_duplicates(
                bin_path=Paths.SVF,
                conf_file_path=run_work_dir / FileNames.CONFIG,
                work_dir_path=run_work_dir,
            )
        else:
            run_work_dir = pathlib.Path(model_find_parameters["run_work_dir"])
        # Get popvoids
        voids, spheres = _pc_postprocessing.get_properties(
            filename=run_work_dir / FileNames.POPFILE
        )
        # Sort voids and spheres ascending in index
        voids.sort_values(by=["id"], inplace=True)
        spheres.sort_values(by=["id"], inplace=True)

        # Get Tracers
        tracers = voids.tracers
        # Get centers of the first sphere
        class0_spheres = spheres[spheres["level"] == 0.0]
        extra = {
            "files_directory_path": run_work_dir,
            "voids": voids,
            "spheres": spheres,
        }
        popcorn_centers = np.array(class0_spheres[["x", "y", "z"]])
        return tracers, popcorn_centers, extra
