#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Rava Jorge Federico, Gualpa Sebastian
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
"""Table reader module."""
import pandas as pd

from . import box


def read_table(path_or_buffer, **kwargs):
    """Input reader.

    Read a table from a file or buffer and returns a `box.Box` object.

    The table must contain 7 columns, with the following format:

    x y z vx vy vz m


    where:

    * `x`, `y`, and `z` are the coordinates of the box center
    * `vx`, `vy`, and `vz` are the velocities of the box center
    * `m` is the mass of the box

    Parameters
    ----------
        path_or_buffer: The path to the file or buffer containing the table.
        **kwargs: Keyword arguments to be passed to `pandas.read_csv()`.

    Returns
    -------
        A `box.Box` object containing the data from the table.

    """
    kwargs.setdefault("sep", r"\s+")
    kwargs.setdefault("usecols", [0, 1, 2, 3, 4, 5, 6])
    data = pd.read_csv(path_or_buffer, **kwargs, header=None)
    col_number = len(data.columns)

    if col_number != 7:
        raise ValueError(
            "There are not enough columns to create the coordinates of a box."
            f"Found {col_number} expected 7"
        )

    check_values = data.notnull().values.all()
    if not check_values:
        raise TypeError(
            f"There are:"
            f"{data.isnull().sum().sum()} null or missing values"
        )

    the_box = box.Box(
        x=data.iloc[:, 0],
        y=data.iloc[:, 1],
        z=data.iloc[:, 2],
        vx=data.iloc[:, 3],
        vy=data.iloc[:, 4],
        vz=data.iloc[:, 5],
        m=data.iloc[:, 6],
    )
    return the_box
