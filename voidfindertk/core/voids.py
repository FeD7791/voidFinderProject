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
"""Void Class."""

# =============================================================================
# IMPORTS
# =============================================================================
import attrs

import numpy as np

from . import center_finder
from . import plot_acc, vsf
from . import radius_finder
from .box import Box
from ..core import cleaners
from ..core import vgcf
from ..utils import Bunch


@attrs.define(frozen=True, repr=False)
class Voids:
    """
    A class to represent and manage voids in a system of particles.

    This class provides methods to analyze voids and the particles within them,
    including retrieving the method used to find the voids, accessing the box
    properties, and finding specific voids associated with tracer particles.

    Parameters
    ----------
    method : str
        The name of the method used to find the voids.
    box : Box
        The Box object containing properties of the box.
    tracers_in_voids : tuple of numpy.ndarray
        A tuple of arrays, where each array contains the IDs of particles
        inside a void.
    centers : array-like
        Coordinates of the centers of the voids.
    extra : dict
        Additional results and information of the run.

    """

    # came from finde and data
    method: str = attrs.field(converter=str)
    box: Box = attrs.field()

    # if end with "_" is calculated
    tracers_in_voids_: tuple = attrs.field(converter=tuple)
    centers_: np.ndarray = attrs.field(converter=np.array)
    extra_: Bunch = attrs.field(converter=lambda e: Bunch("extra", e))

    # plot accessor
    plot: plot_acc.VoidPlotter = attrs.field(
        init=False,
        default=attrs.Factory(plot_acc.VoidPlotter, takes_self=True),
    )

    def __attrs_post_init__(self):
        """Post init method."""
        if len(self.box) <= len(self.tracers_in_voids_):
            raise ValueError(
                "Number of voids can never outnumber number of tracers"
            )

    @property
    def numbers_of_voids_(self):
        """int: Number of voids."""
        return len(self.tracers_in_voids_)

    @property
    def e_(self):
        """dict: Holds extra results and information of the run."""
        return self.extra_

    @property
    def tracers(self):
        """Returns the tracers that belong to a void."""
        return self.tracers_in_voids_

    @property
    def centers(self):
        """Returns the center/s of void/s."""
        return self.centers_

    # REPR ====================================================================
    def __repr__(self):
        """Representation method."""
        return (
            f"<Voids '{self.method}' "
            f"{self.numbers_of_voids_}V, {len(self.box)}T>"
        )

    # utilities ===============================================================

    def filter_by_index(self, index):
        parameters = {
            "method": self.method,
            "box": self.box,
            "tracers_in_voids_": tuple(
                np.array(self.tracers_in_voids_, dtype=object)[index]
            ),
            "centers_": self.centers_[index],
        }

        if "radius" in self.extra_:
            parameters["extra_"] = np.array(self.extra_["radius"])[index]
        else:
            parameters["extra_"] = {}
        return self._create_new_instance(parameters=parameters)

    def void_of(self, tracer):
        """
        Returns indices of voids containing a specific tracer particle.

        Parameters
        ----------
        tracer : int
            ID of the tracer particle to search for in voids.

        Returns
        -------
        numpy.ndarray
            Array of indices of voids containing the tracer particle.

        """
        voids_w_tracer = []
        for idx, void in enumerate(self.tracers_in_voids_):
            if tracer in void:
                voids_w_tracer.append(idx)
        return np.array(voids_w_tracer)

    def effective_radius(self, method="density", **kwargs):
        """
        Computes the effective radius of voids using the given method.

        Parameters
        ----------
        method : {'density', 'extra', 'volume'}, default='density'
            The method used to calculate the effective radius.
        **kwargs : keyword arguments
            Additional keyword arguments passed to the method used to calculate
            the effective radius.

        Returns
        -------
        numpy.ndarray
            An array containing the calculated effective radii of the voids.
        """
        rad_method = radius_finder.get_radius_searcher(method=method)
        return rad_method(
            centers=self.centers_, box=self.box, extra=self.extra_, **kwargs
        )

    def void_size_function(
        self,
        *,
        radius,
        delta=-0.8,
        n_step1=2,
        n_step2=10,
        scale_1_num_samples=7,
        scale_2_num_samples=2,
    ):
        """Calculate the void size distribution based on the effective radius.

        This function computes the log of radius, count of voids, and delta
        values using the provided parameters for scaling and sampling.

        Parameters
        ----------
        n_step1 : int, optional
            Steps for the first scale. Default is 2. The lower the more number
            of tracers (See vsf.void_size_function)
        n_step2 : int, optional
            Steps for the second scale. Default is 10.
        scale_1_num_samples : int, optional
            Number of the first part of bins. Default is 7.
        scale_2_num_samples : int, optional
            Numbers of the second part of bins. Default is 2.
        **kwargs : keyword arguments
            Additional keyword arguments passed to the `effective_radius`
            method.

        Returns
        -------
        tuple
            A tuple containing:
            - log_of_radius : array
                The logarithm of the effective radius values.
            - count : array
                The count of voids corresponding to each radius.
            - delta : array
                The delta values for the void size distribution.
        """
        log_of_radius, count, delta = vsf.void_size_function(
            radius=radius,
            delta=delta,
            box=self.box,
            n_step1=n_step1,
            n_step2=n_step2,
            scale_1_num_samples=scale_1_num_samples,
            scale_2_num_samples=scale_2_num_samples,
        )

        return log_of_radius, count, delta

    vsf = void_size_function

    def _cleaner(
        self,
        cleaner_method="overlap",
        rad_min_max=[0.0, 100.0],
        **kwargs,
    ):
        """
        Cleans voids based on a cleaner method and given radius.

        This method uses the specified cleaner method to clean the voids and
        applies a filtering process based on the given radius.

        Parameters
        ----------
        centers : array-like
            The coordinates of the centers of the voids.
        radius : array-like
            The radii of the voids.
        cleaner_method : {'overlap', 'cbl'}, default='overlap'
            The method used to clean the voids.
        rad_min_max : list, optional
            A list containing the minimum and maximum radius for filtering
            voids. Default is [0.0, 100.0].
        **kwargs : keyword arguments
            Additional keyword arguments passed to the cleaner method.

        Returns
        -------
        Voids
            A new instance of the Voids class with cleaned voids.
        """

        cleaner = cleaners.get_cleaner(cleaner_method=cleaner_method)
        centers_, radius_ = cleaner(**kwargs)
        # Post cleaning
        idx = np.where(
            (radius_ < rad_min_max[1]) & (radius_ > rad_min_max[0])
        )[0]
        # Build new void
        parameters = {
            "method": self.method,
            "box": self.box,
            "tracers_in_voids_": tuple(
                [self.tracers_in_voids_[i] for i in idx]
            ),
            "centers_": centers_[idx],
            "extra_": {"radius": radius_[idx]},
        }
        return self._create_new_instance(parameters=parameters)

    @classmethod
    def _create_new_instance(cls, parameters):
        """
        Creates a new instance of the Voids class with the given parameters.

        Parameters
        ----------
        parameters : dict
            A dictionary containing the parameters to initialize the new Voids
            object.

        Returns
        -------
        Voids
            A new instance of the Voids class initialized with the given
            parameters.
        """
        return cls(
            method=parameters["method"],
            box=parameters["box"],
            tracers_in_voids_=parameters["tracers_in_voids_"],
            centers_=parameters["centers_"],
            extra_=parameters["extra_"],
        )

    def find_radius_and_clean(
        self, cleaner_method="overlap", radius_method="density", **kwargs
    ):
        """
        Finds the radius of voids and cleans them using the specified methods.

        This method calculates the effective radius of voids using the
        specified radius method, applies a cleaning process, and returns a new
        instance of the Voids class with cleaned voids.

        Parameters
        ----------
        cleaner_method : str, {'overlap','cbl'}, default='overlap'
            The cleaner method to use.
        radius_method : str, {'density','extra','volume'}
            The method to use for calculating the effective radius. Default is
            "default".
        **kwargs : keyword arguments
            Additional keyword arguments passed to the effective radius and
            cleaner methods.

        Returns
        -------
        Voids
            A new instance of the Voids class with cleaned voids and calculated
            radii.
        """

        # 1) Find radius by method radius_method
        radius = self.effective_radius(method=radius_method, **kwargs)

        # 2) Clean found radius
        # 2.1) Pre cleaning
        idx = np.where(radius > 0)[0]
        radius = radius[idx]
        centers = self.centers_[idx]

        # Set default params
        kwargs.setdefault("box", self.box)
        kwargs.setdefault("extra", self.extra_)
        kwargs["center"] = centers
        kwargs["radius"] = radius

        # 2.2) Use class cleaner method
        new_void = self._cleaner(cleaner_method=cleaner_method, **kwargs)
        return new_void

    def find_centers(
        self,
        n_neighbors=3,
        threshold=0.8,
        n_tracers_threshold=[0, 1000],
        n_jobs=1,
        batch_size=10,
    ):
        """
        Finds and calculates the centers of voids.

        This method computes the centers of the voids based on a set of tracers
        that are inside voids.

        Returns
        -------
        Voids
            A new instance of the Voids class with the new centers, the
            tracers inside each void are the same.
        """
        # Get number of tracers in each void
        n_tracers = np.array(list(map(len, self.tracers_in_voids_)))

        # Filter centers by number of tracers
        idx = np.where(
            (n_tracers > n_tracers_threshold[0])
            & (n_tracers < n_tracers_threshold[1])
        )[0]

        # Briefly map tracers to numpy array to perform filter
        _tracers = tuple([self.tracers_in_voids_[i] for i in idx])

        new_centers = center_finder.center_calculator(
            box=self.box,
            tracers_in_voids=_tracers,
            n_neighbors=n_neighbors,
            threshold=threshold,
            n_jobs=n_jobs,
            batch_size=batch_size,
        )
        # Build new void
        parameters = {
            "method": self.method,
            "box": self.box,
            "tracers_in_voids_": _tracers,
            "centers_": new_centers,
            "extra_": {},
        }
        return self._create_new_instance(parameters=parameters)

    def void_galaxy_corr(self, max_radius_search=30, delta_rad=1, n_jobs=1):
        """
        Calculates the void-galaxy cross correlation function.

        The search is performed by calculating the tracers inside of concentric
        spherical shells with width [rad ,rad+delta_rad] up to
        max_radius_search distance around void centers.

        Parameters
        ----------
        max_radius_search : float
            Maximum radius of the search.
        delta_rad : float, default=1
            Width in r of the concentric shells.
        n_jobs : int , default=1
            Measure of the number of cores used to perform a parallel search.
            Use -1 to exhaust the cores.

        Returns
        -------
        vgcf : list
            Values of the vgcf.
        """
        return vgcf.vgcf_statistic(
            centers=self.centers_,
            box=self.box,
            max_rad=max_radius_search,
            delta_r=delta_rad,
            n_jobs=n_jobs,
        )
