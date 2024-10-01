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
"""DIVE: Detection, Identificacion and Voids Examination.

DIVE is a reimplementation of VIDE void detection algoritm in pure Python on
top of the classic ZOBOV code.
"""
# =============================================================================
# IMPORTS
# =============================================================================
from ._dive import DiveVF

# =============================================================================
# ALL
# =============================================================================
__all__ = ["DiveVF"]
