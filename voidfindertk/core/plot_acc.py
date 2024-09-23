#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Rava Jorge Federico, Gualpa Sebastian
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.

# =============================================================================
# DOCS
# =============================================================================
"""
Module for plotting data related to boxes and voids using Matplotlib and
Seaborn.

This module contains classes that facilitate the visualization of data 
from box and void objects through various plotting functions.

Classes
-------
BoxPlotter
    A class for plotting histograms of data contained in a box.

VoidPlotter
    A class for plotting void size functions.

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
        self._box = box

    def hist2d(self, x, y, *, ax=None, **kwargs):
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
        self._voids = voids

    def __getattr__(self, kind):
        voids = self._voids
        return getattr(voids.box.plot, kind)

    def __setstate__(self, state):
        self.__dict__.update(state)

    def void_size_function(self, *, ax=None, vsf_kws=None, **kwargs):
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
