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

"""Test for voidfindertk.utils.bunch."""


# =============================================================================
# IMPORTS
# =============================================================================

import copy
import pickle

import pytest

from voidfindertk.utils import bunch


# =============================================================================
# TEST Bunch
# =============================================================================


def test_bunch_creation():
    md = bunch.Bunch("foo", {"alfa": 1})
    assert md["alfa"] == md.alfa == 1
    assert len(md) == 1


def test_bunch_creation_empty():
    md = bunch.Bunch("foo", {})
    assert len(md) == 0


def test_bunch_key_notfound():
    md = bunch.Bunch("foo", {"alfa": 1})
    assert md["alfa"] == md.alfa == 1
    with pytest.raises(KeyError):
        md["bravo"]


def test_bunch_attribute_notfound():
    md = bunch.Bunch("foo", {"alfa": 1})
    assert md["alfa"] == md.alfa == 1
    with pytest.raises(AttributeError):
        md.bravo


def test_bunch_iter():
    md = bunch.Bunch("foo", {"alfa": 1})
    assert list(iter(md)) == ["alfa"]


def test_bunch_repr():
    md = bunch.Bunch("foo", {"alfa": 1})
    assert repr(md) == "<foo {'alfa'}>"


def test_bunch_dir():
    md = bunch.Bunch("foo", {"alfa": 1})
    assert "alfa" in dir(md)


def test_bunch_deepcopy():
    md = bunch.Bunch("foo", {"alfa": 1})
    md_c = copy.deepcopy(md)

    assert md is not md_c
    assert md._name == md_c._name  # string are inmutable never deep copy
    assert md._data == md_c._data and md._data is not md_c._data


def test_bunch_copy():
    md = bunch.Bunch("foo", {"alfa": 1})
    md_c = copy.copy(md)

    assert md is not md_c
    assert md._name == md_c._name
    assert md._data == md_c._data and md._data is md_c._data


def test_bunch_setstate():
    md = bunch.Bunch("foo", {"alfa": 1})
    reloaded = pickle.loads(pickle.dumps(md))
    assert md == reloaded


def test_bunch_todict():
    md = bunch.Bunch("foo", {"alfa": 1})
    assert md.to_dict() == {"alfa": 1}
