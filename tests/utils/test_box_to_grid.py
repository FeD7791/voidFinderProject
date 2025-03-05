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

"""Test for bot to grispy grid parsing."""

# =============================================================================
# IMPORTS
# =============================================================================

import grispy as gsp

from voidfindertk.utils import box_to_grid

# =============================================================================
# TESTS
# =============================================================================


def test_get_grispy_grid_from_box_default(mkbox):
    """Test the default behavior of get_grispy_grid_from_box."""

    mock_box = mkbox(seed=42, size=1000)
    # Call the function with the mock box and no additional arguments
    grid = box_to_grid.get_grispy_grid_from_box(mock_box)

    # Assert the grid is of type gsp.GriSPy
    assert isinstance(grid, gsp.GriSPy)

    # Check if the grid's periodicity is correctly set for all axes (x, y, z)
    periodic = {
        0: (mock_box.min_, mock_box.max_),
        1: (mock_box.min_, mock_box.max_),
        2: (mock_box.min_, mock_box.max_),
    }
    grid.set_periodicity(periodic, inplace=True)

    # Verify if default grid parameters are applied correctly
    assert grid.N_cells == 64  # Default value for N_cells
    assert grid.copy_data is False  # Default value for copy_data


def test_get_grispy_grid_from_box_with_custom_params(mkbox):

    mock_box = mkbox(seed=42, size=1000)
    """Test the behavior of get_grispy_grid_from_box with custom parameters."""

    # Custom parameters to pass to the function
    grispy_kwargs = {"N_cells": 128, "copy_data": True}

    # Call the function with the mock box and custom arguments
    grid = box_to_grid.get_grispy_grid_from_box(mock_box, **grispy_kwargs)

    # Assert the grid is of type gsp.GriSPy
    assert isinstance(grid, gsp.GriSPy)

    # Check if the grid's periodicity is correctly set for all axes (x, y, z)
    periodic = {
        0: (mock_box.min_, mock_box.max_),
        1: (mock_box.min_, mock_box.max_),
        2: (mock_box.min_, mock_box.max_),
    }
    grid.set_periodicity(periodic, inplace=True)

    # Verify if custom grid parameters are applied correctly
    assert grid.N_cells == 128  # Custom value for N_cells
    assert grid.copy_data is True  # Custom value for copy_data
