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

""" test for voidfindertk.utils.context_managers

"""

# =============================================================================
# IMPORTS
# =============================================================================

import contextlib
import os
import tempfile


import pytest

# =============================================================================
# FUNCTIONS
# =============================================================================


@pytest.fixture
def temp_dir():
    """Fixture to create a temporary directory for testing."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


def test_chdir_with_temp_directory(temp_dir):
    # Get the original current working directory
    original_cwd = os.getcwd()

    # Use the context manager to change the directory
    with contextlib.chdir(temp_dir):
        # Assert that the current directory is now the temp directory
        assert os.getcwd() == temp_dir

    # After exiting the context manager, assert that the current directory has
    # returned to the original
    assert os.getcwd() == original_cwd
