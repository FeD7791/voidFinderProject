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

"""C objects builder."""

# =============================================================================
# IMPORTS
# =============================================================================

import os
import subprocess

from setuptools import setup
from setuptools.command.build_py import build_py

# =============================================================================
# BUILDER
# =============================================================================


class BuildWithMake(build_py):
    """

    Custom setuptools command to compile C libraries using a Makefile before \
    building Python packages.

    This class overrides the default `build_py` command to:
    1. Run `make` in the specified subdirectory (`voidfindertk/zobov/`).
    2. Ensure the compiled `.so` files remain in their original location.
    3. Proceed with the standard Python build process.

    Notes
    -----
    - The Makefile must generate shared libraries (`.so` files) in
    `voidfindertk/zobov/`.
    - The `.so` files are automatically included in the package via
    `pyproject.toml`.

    Examples
    --------
    When the package is installed via `pip install .`, this command:
    1. Executes `make` in `voidfindertk/zobov/`.
    2. Builds the Python package while preserving the `.so` files in `zobov/`.

    """

    def run(self):
        """Runs make in voidfindertk/zobov/."""
        make_dir = os.path.join("voidfindertk", "zobov")
        subprocess.check_call(["make"], cwd=make_dir)
        super().run()  # Proceed with normal Python build


# Minimal setup() - metadata is in pyproject.toml
setup(cmdclass={"build_py": BuildWithMake})
