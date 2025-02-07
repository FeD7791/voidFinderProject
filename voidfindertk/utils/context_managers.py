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

"""Change dir context manager."""

# =============================================================================
# IMPORTS
# =============================================================================


import contextlib
import os

# =============================================================================
# FUNCTIONS
# =============================================================================

try:
    chdir = contextlib.chdir
except AttributeError:

    class chdir(contextlib.AbstractContextManager):
        """Non thread-safe context manager to change the current working \
        directory."""

        def __init__(self, path):
            self.path = path
            self._old_cwd = []

        def __enter__(self):
            """Context manager enter method."""
            self._old_cwd.append(os.getcwd())
            os.chdir(self.path)

        def __exit__(self, *excinfo):
            """Context manager exit method."""
            os.chdir(self._old_cwd.pop())
