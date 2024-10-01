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
import attr

import numpy as np

from . import _postprocessing, _wrapper
from ..svf_popcorn import FileNames, Paths, SVFPopCorn


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

    _shot_noise_threshold = attr.field(default=20)

    @property
    def shot_noise_threshold(self):
        """Return the shot noise threshold value."""
        return self._shot_noise_threshold

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
        tracers_in_voids, centers, extra = super().build_voids(
            model_find_parameters
        )
        run_work_dir = extra["files_directory_path"]

        # Before continuing minradius must be re-configured
        _wrapper.read_and_modify_config(
            config_file_path=run_work_dir / FileNames.CONFIG,
            section="INPUT_PARAMS",
            parameter="MINRADIUS",
            new_value=str(self._shot_noise_threshold),
        )
        _wrapper.popcorn_void_finder(
            mpi_flags=self._mpi_flags,
            bin_path=Paths.SVF,
            conf_file_path=run_work_dir / FileNames.CONFIG,
            work_dir_path=run_work_dir,
        )
        _wrapper.compute_intersects(
            bin_path=Paths.SVF,
            conf_file_path=run_work_dir / FileNames.CONFIG,
            work_dir_path=run_work_dir,
        )
        _wrapper.clean_duplicates(
            bin_path=Paths.SVF,
            conf_file_path=run_work_dir / FileNames.CONFIG,
            work_dir_path=run_work_dir,
        )
        # Get popvoids
        popcorn_void_properties = _postprocessing.read_pop(
            filename=run_work_dir / FileNames.POPFILE
        )
        # Process output
        tracers = []
        # Get centers of the first sphere
        # By now, the centers are the center of the spheres of heriarchy level
        # zero
        x = np.array(
            [e["x"][0] for e in popcorn_void_properties["pop"]],
            dtype=np.float32,
        )
        y = np.array(
            [e["y"][0] for e in popcorn_void_properties["pop"]],
            dtype=np.float32,
        )
        z = np.array(
            [e["z"][0] for e in popcorn_void_properties["pop"]],
            dtype=np.float32,
        )
        centers = np.stack([x, y, z], axis=1)
        extra = {
            "files_directory_path": run_work_dir,
            "r_eff": np.array([e for e in popcorn_void_properties["reff"]]),
        }
        return tracers, centers, extra
