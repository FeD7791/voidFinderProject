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

"""test for voidfindertk.utils.accabc

"""


# =============================================================================
# IMPORTS
# =============================================================================


import numpy as np

import pytest

from voidfindertk.utils.accabc import AccessorABC


# =============================================================================
# TEST CLASSES
# =============================================================================


def test_AccessorABC():
    class FooAccessor(AccessorABC):
        _default_kind = "zaraza"

        def __init__(self, v):
            self._v = v

        def zaraza(self):
            return self._v

    acc = FooAccessor(np.random.random())
    assert acc("zaraza") == acc.zaraza() == acc()


def test_AccessorABC_no__default_kind():
    with pytest.raises(TypeError):

        class FooAccessor(AccessorABC):
            pass

    with pytest.raises(TypeError):
        AccessorABC()


def test_AccessorABC_invalid_kind():
    class FooAccessor(AccessorABC):
        _default_kind = "zaraza"

        def __init__(self):
            self.dont_work = None

        def _zaraza(self):
            pass

    acc = FooAccessor()

    with pytest.raises(ValueError):
        acc("_zaraza")

    with pytest.raises(ValueError):
        acc("dont_work")
