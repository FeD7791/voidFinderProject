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

"""Test for voidfindertk.utils.make_workdir."""

# =============================================================================
# IMPORTS
# =============================================================================

from unittest import mock


from voidfindertk.utils import make_workdir

# =============================================================================
# TESTS
# =============================================================================


def test_make_workdir():

    with mock.patch("tempfile.mkdtemp", return_value="path"):
        path = make_workdir.create_run_work_dir(workdir_path="wdpath")

    assert str(path) == "path"
