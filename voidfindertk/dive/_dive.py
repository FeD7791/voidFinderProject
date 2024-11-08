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
"""Module that holds functions and methods that are used to run the DIVE \
void finder."""
# =============================================================================
# IMPORTS
# =============================================================================
from . import _postprocessing
from ..zobov import Names, ZobovVF


class _Files:
    """
    A class containing file name constants for output files.

    This class defines several constants that specify the names of
    various output files used in the processing of void and tracer data.

    Attributes
    ----------
    VOL_FILE_RAW : str
        The name of the raw volume file, formatted with the output name.
    ADJ_FILE_RAW : str
        The name of the raw adjacency file, formatted with the output name.
    XYZ_R_EFF_VOIDS_FILE : str
        The name of the file containing void coordinates and effective radius.
    XYZ_TRACERS_FILE : str
        The name of the file containing tracer coordinates.
    CLEANED_CATALOGUE : str
        The name of the cleaned catalogue file containing processed void data.

    Notes
    -----
    The `Names.OUTPUT_VOZINIT` should be defined elsewhere in the code to
    provide a valid output name used in the volume and adjacency file names.
    """

    VOL_FILE_RAW = f"vol{Names.OUTPUT_VOZINIT}.dat"
    ADJ_FILE_RAW = f"adj{Names.OUTPUT_VOZINIT}.dat"
    XYZ_R_EFF_VOIDS_FILE = "xyz_r_eff_file.txt"
    XYZ_TRACERS_FILE = "xyz_tracers.txt"
    CLEANED_CATALOGUE = "cleaned_catalogue.txt"


class DiveVF(ZobovVF):
    """
    Class for processing voids using the ZOBOV Finder algorithm.

    This class extends the `ZobovVF` class to implement specific parameters
    and methods for building voids from tracer data.

    Parameters
    ----------
    ratio : float, optional
        The ratio used in void cleaning (default is 0.1).
    initial_radius : bool, optional
        Flag to indicate if the initial radius should be used (default is
        True).
    delta_r : list of float, optional
        Range of delta radii for processing (default is [17.0, 150.0]).
    threshold : float, optional
        Threshold value for void detection (default is 0.3).
    overlap_criterion : bool, optional
        Flag to determine if overlap criterion is applied (default is True).
    **kwargs :
        Additional keyword arguments passed to the parent class.


    """

    def __init__(
        self,
        *,
        ratio=0.1,
        initial_radius=True,
        delta_r=[17.0, 150.0],
        threshold=0.3,
        overlap_criterion=True,
        **kwargs,
    ):
        """Init class method."""
        super().__init__(**kwargs)
        self._ratio = ratio
        self._initial_radius = initial_radius
        self._delta_r = delta_r
        self._threshold = threshold
        self._overlap_criterion = overlap_criterion

    @property
    def ratio(self):
        """Get the ratio used in void cleaning."""
        return self._ratio

    @property
    def initial_radius(self):
        """Get initial radius."""
        return self._initial_radius

    @property
    def delta_r(self):
        """Get the range of delta radii for processing."""
        return self._delta_r

    @property
    def threshold(self):
        """Get the threshold for void detection."""
        return self._threshold

    @property
    def overlap_criterion(self):
        """Get overlap criterion is applied."""
        return self._overlap_criterion

    def build_voids(self, model_find_parameters):
        """
        Method that postprocess the outputs from ZOBOV Finder.

        Parameters
        ----------
            model_find_parameters : dict
            Dict object with the following keys:
                run_work_dir : pathlib.Path
                path to run work directory inside workdir.
                box : object
                box object
        """
        tracers_in_voids, zobov_centers, extra = super().build_voids(
            model_find_parameters
        )
        # Get working directory
        run_work_dir = extra["files_directory_path"]
        # Get void properties
        void_properties = extra["void_properties"]
        # Get box
        box = model_find_parameters["box"]
        # Get tracer volumes
        tracer_volumes = _postprocessing.read_volume_file(
            filename=run_work_dir / _Files.VOL_FILE_RAW
        )
        # get center and radii
        radii, centers = _postprocessing.get_center_and_radii(
            void_properties=void_properties,
            tracer_volumes=tracer_volumes,
            tracers_in_voids=tracers_in_voids,
            box=box,
        )
        # Save to file
        # 1) center and radii
        _postprocessing.save_r_eff_center(
            centers=centers,
            r_eff=radii,
            path=str(run_work_dir / _Files.XYZ_R_EFF_VOIDS_FILE),
        )
        # 2) xyz box
        _postprocessing.save_xyz_tracers(
            box=box, path=str(run_work_dir / _Files.XYZ_TRACERS_FILE)
        )
        # Perform cleaning
        _postprocessing.cbl_cleaner(
            file_voids=str(run_work_dir / _Files.XYZ_R_EFF_VOIDS_FILE),
            file_tracers=str(run_work_dir / _Files.XYZ_TRACERS_FILE),
            ratio=self._ratio,
            initial_radius=self._initial_radius,
            delta_r=self._delta_r,
            threshold=self._threshold,
            output_path=str(run_work_dir / _Files.CLEANED_CATALOGUE),
            ol_crit=self._overlap_criterion,
        )
        # Get tracers in voids and center
        clean_tracers, clean_centers = _postprocessing.get_tracers_and_center(
            box=box,
            cbl_cleaned_path=str(run_work_dir / _Files.CLEANED_CATALOGUE),
        )
        # Updating extra with void_properties
        extra["void_properties"] = _postprocessing.get_dive_void_properties(
            cleaned_catalogue_path=run_work_dir / _Files.CLEANED_CATALOGUE
        )
        return tuple(clean_tracers), clean_centers, extra
