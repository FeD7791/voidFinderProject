#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Rava Jorge Federico, Gualpa Sebastian
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
""" ABC for Void Search """
from abc import ABC, abstractmethod

from .voids import Voids


class VoidFinderABC(ABC):
    """
    Abstract base class for finding voids in a given box.

    This class defines the interface for subclasses to implement the
    functionality to preprocess data, model the findings, and build voids.

    Methods
    -------
    find(box)
        Main method to find voids in the provided box.

    preprocess(box)
        Abstract method to preprocess the input box. Must be implemented by
        subclass.

    model_find(preprocess_parameters)
        Abstract method to model the findings from the preprocessed parameters.
        Must be implemented by subclass.

    build_voids(model_find_parameters)
        Abstract method to build voids from the model find parameters.
        Must be implemented by subclass.
    """
    def __init__(self):
        """Initialize the VoidFinderABC instance."""
        pass

    def find(self, box):
        """
        Find voids in the provided box.

        Parameters
        ----------
        box : object
            The input data structure representing the box of tracers in which
            voids are to be found.

        Returns
        -------
        Voids
            An instance of the Voids class containing the voids information,
            including tracers in voids, their centers, and any extra
            information.
        """
        preprocess_parameters = self.preprocess(box)
        model_find_parameters = self.model_find(preprocess_parameters)
        tracers_in_voids, centers, extra = self.build_voids(
            model_find_parameters
        )

        voids = Voids(
            method=type(self).__name__,
            box=box,
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
            The parameters resulting from preprocessing the input box.
        """
        pass

    @abstractmethod
    def model_find(self, preprocess_parameters):
        """
        Execute the search for voids in the Box.

        Parameters
        ----------
        preprocess_parameters : object
            The parameters resulting from the preprocessing step.

        Returns
        -------
        model_find_parameters : object
            The parameters resulting from the model finding process.
        """
        pass

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
        pass
