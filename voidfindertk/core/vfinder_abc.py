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
"""ABC for Void Search."""

# =============================================================================
# IMPORTS
# =============================================================================
from abc import ABC, abstractmethod

from .voids import Voids


class VoidFinderABC(ABC):
    """
    Abstract base class for finding voids in a given box.

    This class defines the interface for subclasses to implement the
    functionality to preprocess data, model the findings, and build voids.

    """

    def __init__(self):
        """Initialize the VoidFinderABC instance."""
        pass

    def find(self, box, box_copy=False):
        """
        Find voids in the provided box.

        Parameters
        ----------
        box : object
            The input data structure representing the box of tracers in which
            voids are to be found.
        box_copy : bool, optional
            Whether to copy the box before preprocessing it. Default is False.

        Returns
        -------
        Voids
            An instance of the Voids class containing the voids information,
            including tracers in voids, their centers, and any extra
            information.
        """
        procesed_box = self.preprocess(box)
        model_find_parameters = self.model_find(procesed_box)
        tracers_in_voids, centers, extra = self.build_voids(
            model_find_parameters
        )

        voids = Voids(
            method=type(self).__name__,
            box=procesed_box.copy() if box_copy else procesed_box,
            tracers_in_voids_=tracers_in_voids,
            centers_=centers,
            extra_=extra,
        )

        return voids

    @abstractmethod
    def preprocess(self, box):
        """
        Preprocess the input box if needed.

        Parameters
        ----------
        box : object
            The input data structure to preprocess.

        Returns
        -------
        preprocess_parameters : object
            Returns Box object with procesed parameters.
        """
        pass  # pragma: no cover

    @abstractmethod
    def model_find(self, procesed_box):
        """
        Execute the search for voids in the Box.

        Parameters
        ----------
        procesed_box : object
            Box object resulting from the preprocessing step.

        Returns
        -------
        model_find_parameters : object
            The parameters resulting from the model finding process.
        """
        pass  # pragma: no cover

    @abstractmethod
    def build_voids(self, model_find_parameters):
        """
        Build voids from the model find parameters.

        Parameters
        ----------
        model_find_parameters : object
            The parameters resulting from the model finding process.

        Returns
        -------
        tracers_in_voids : object
            Information about tracers within the identified voids.
        centers : object
            The centers of the identified voids.
        extra : object
            Any additional information generated during void construction.
        """
        pass  # pragma: no cover
