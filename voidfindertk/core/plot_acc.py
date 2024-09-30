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
"""
Module for plotting data related to boxes and voids.

This module contains classes that facilitate the visualization of data
from box and void objects through various plotting functions. The module uses
Matplotlib and Seaborn.


"""

# =============================================================================
# IMPORTS
# =============================================================================

import matplotlib.pyplot as plt

import seaborn as sns

from ..utils import accabc



class BoxPlotter(accabc.AccessorABC):
    """Plotter: Plotter object that plots the boxes."""

    _default_kind = "hist2d"

    def __init__(self, box):
        """
        Initialize the BoxPlotter object for plotting boxes.

        This constructor initializes a BoxPlotter instance that can create
        2D histogram plots based on the data contained within the specified
        box.

        Parameters
        ----------
        box : object
            An object that contains the data to be plotted. The object is
            expected to have attributes corresponding to the data series
            for the x and y axes used in the plots.

        Notes
        -----
        The `box` parameter should have attributes accessible via the dot
        notation that represent the data to be visualized. This is crucial
        for the proper functioning of the `hist2d` method.
        """
        self._box = box

    def hist2d(self, x, y, *, ax=None, **kwargs):
        """
        Create a 2D histogram plot.

        This method generates a 2D histogram plot of the specified x and
        y data from the box. It uses seaborn's histplot function for
        visualization.

        Parameters
        ----------
        x : str
            The name of the attribute in the box for the x-axis data.
        y : str
            The name of the attribute in the box for the y-axis data.
        ax : matplotlib.axes.Axes, optional
            The axes on which to plot. If None, the current axes are used.
        **kwargs :
            Additional keyword arguments passed to the seaborn histplot.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The axes with the histogram plot.

        Notes
        -----
        The x and y attributes in the box should have a `unit` attribute
        for labeling the axes correctly.

        """
        ax = plt.gca() if ax is None else ax

        xvalues = getattr(self._box, x)
        yvalues = getattr(self._box, y)

        sns.histplot(x=xvalues, y=yvalues, ax=ax, **kwargs)

        ax.set_xlabel(f"{x} ({xvalues.unit})")
        ax.set_ylabel(f"{y} ({yvalues.unit})")

        return ax


class VoidPlotter(accabc.AccessorABC):
    """Plotter: Plotter object that plots the voids."""

    _default_kind = "void_size_function"

    def __init__(self, voids):
        """
        Initialize the VoidPlotter object for plotting voids.

        This constructor initializes a VoidPlotter instance that can create
        visualizations related to voids, specifically the void size function.

        Parameters
        ----------
        voids : object
            An object containing void data to be plotted. It is expected
            to have a `box` attribute with a `plot` method for accessing
            plotting functionalities.

        Notes
        -----
        The `voids` parameter should have attributes that allow for
        plotting operations. This is crucial for the proper functioning
        of the plotting methods in this class.
        """
        self._voids = voids

    def __getattr__(self, kind):
        """Get attribute access to kind."""
        voids = self._voids
        return getattr(voids.box.plot, kind)

    def __setstate__(self, state):
        """Set state."""
        self.__dict__.update(state)

    def void_size_function(self, *, ax=None, vsf_kws=None, **kwargs):
        """
        Create a plot of the void size function.

        This method generates a plot for the void size function using data
        derived from the voids. It visualizes the relationship between the
        logarithm of the radius and the density of voids.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The axes on which to plot. If None, the current axes are used.
        vsf_kws : dict, optional
            Additional keyword arguments passed to the void size function
            calculation.
        **kwargs :
            Additional keyword arguments passed to the seaborn lineplot.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The axes with the void size function plot.

        Notes
        -----
        The x-axis is logarithmically scaled, and the y-axis represents
        the density of voids. The plot includes a grid and is titled with
        the density contrast.
        """
        # confs
        vsf_kws = {} if vsf_kws is None else vsf_kws

        # if no axis get the default
        ax = plt.gca() if ax is None else ax

        # get the vsf
        x, y, delta = self._voids.vsf(**vsf_kws)

        kwargs.setdefault("marker", "o")
        kwargs.setdefault("label", "Void Size Function")
        sns.lineplot(x=x, y=y, ax=ax, **kwargs)

        ax.set_xlabel(r"$log_{10}(R)$" f"{x.unit}")
        ax.set_ylabel(r"$\frac{1}{V} \frac{dN_v}{dlnR_v}$")

        ax.set_yscale("log")

        ax.set_title(f"Void Size Function\n Density Contrast {delta}")

        ax.grid(True)

        return ax

    vsf = void_size_function

    def void_over_hist2d(self, x, y, *, ax=None):
        ax = self.hist2d(x=x, y=y, ax=ax)
        ax.axvline(500, color="k")
