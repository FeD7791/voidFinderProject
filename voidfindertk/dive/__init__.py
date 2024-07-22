#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the
#   Void-Finder-Toolkit Project (https://github.com/FeD7791/voidFinderProject).
# Copyright (c) 2023-2024, Bustillos, J.F.; Cabral, Juan
# License: BSD 3-Clause
# Full Text:
#   https://raw.githubusercontent.com/FeD7791/voidFinderProject/dev/LICENSE.txt

"""DIVE: Detection, Identificacion and Voids Examination

DIVE is a reimplementation of VIDE void detection algoritm in pure Python on
top of the classic ZOBOV code.

"""

from ._processing import get_center_and_radii

__all__ = ["get_center_and_radii"]