###############################################################################
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Federico, Gualpa Sebastian, Cabral Juan
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
###############################################################################
import pytest

from voidfindertk import box
from voidfindertk.io import read_table


def test_read_table_type_output(random_buffer):
    buffer = random_buffer()
    assert isinstance(read_table(buffer), box.Box)


def test_nullvalues_on_read(random_buffer):
    buffer = random_buffer(empty_row=True)
    with pytest.raises(TypeError):
        read_table(buffer)


def test_number_of_columns(random_buffer):
    buffer = random_buffer()
    with pytest.raises(ValueError):
        read_table(buffer, usecols=[1, 2, 3])
