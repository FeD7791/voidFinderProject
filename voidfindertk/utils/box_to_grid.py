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

"""Parses Box object to Grispy Grid."""

# =============================================================================
# IMPORTS
# =============================================================================

import grispy as gsp

import numpy as np

# =============================================================================
# FUNCTIONS
# =============================================================================


def get_grispy_grid_from_box(box, **grispy_kwargs):
    """Gets a periodical grispy grid from box."""
    grispy_kwargs.setdefault("N_cells", 64)
    grispy_kwargs.setdefault("copy_data", False)

    # Build a grispy grid using the box
    x = box.arr_.x
    y = box.arr_.y
    z = box.arr_.z
    xyz = np.column_stack((x, y, z))
    grid = gsp.GriSPy(
        xyz,
        copy_data=grispy_kwargs["copy_data"],
        N_cells=grispy_kwargs["N_cells"],
    )
    periodic = {
        0: (box.min(), box.max()),
        1: (box.min(), box.max()),
        2: (box.min(), box.max()),
    }
    grid.set_periodicity(periodic, inplace=True)
    return grid
