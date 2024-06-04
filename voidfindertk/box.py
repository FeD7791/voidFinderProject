#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Rava Jorge Federico, Gualpa Sebastian
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
"""Class Box object constructor."""
import math

from astropy import units as u

import numpy as np

import uttr


@uttr.s(repr=False, frozen=True, cmp=False)
class Box:
    """Box Class.

    Class used to describe a set of points (x,y,z) alongside with its
    velocities (vx,vy,vz)

    Attributes
    ----------
    x : numpy.ndarray
    y : numpy.ndarray
    z : numpy.ndarray
        (x,y,z) array of position elements
    vx : numpy.ndarray
    vy : numpy.ndarray
    vz : numpy.ndarray
        (vx,vy,vz) array of velocity elements

    Methods
    -------
    build_df(self)
        Transforms the set into a pandas dataframe

    _check_length(self)
        Verifies that the lenght of the inputs are the same

    slice(self,n,parameter)
        Calls the function Slicer to slice the set in n parts
    """

    x = uttr.ib(converter=np.array, unit=u.Mpc)
    y = uttr.ib(converter=np.array, unit=u.Mpc)
    z = uttr.ib(converter=np.array, unit=u.Mpc)
    vx = uttr.ib(converter=np.array, unit=u.Mpc / u.second)
    vy = uttr.ib(converter=np.array, unit=u.Mpc / u.second)
    vz = uttr.ib(converter=np.array, unit=u.Mpc / u.second)
    m = uttr.ib(converter=np.array, unit=u.M_sun)

    _len = uttr.ib(init=False)

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        lengths = set(())
        for e in (self.x, self.y, self.z, self.vx, self.vy, self.vz):
            lengths.add(len(e))

        if len(lengths) != 1:
            raise ValueError("Arrays should be of the same size")

        super().__setattr__("_len", lengths.pop())

        # check if the box is a cube
        box_side = set(
            (
                math.ceil(np.max(self.x.value)),
                math.ceil(np.max(self.y.value)),
                math.ceil(np.max(self.z.value)),
            )
        )
        if len(box_side) != 1:
            raise ValueError(
                "Not a cube:"
                + f" xmax: {math.ceil(np.max(self.x.value))}"
                + f" ymax: {math.ceil(np.max(self.y.value))}"
                + f" zmax: {math.ceil(np.max(self.z.value))}"
            )

    def __len__(self):
        """Length method.

        Returns
        -------
            int
                the number of elements in the box
        """
        return self._len

    def __eq__(self, other):
        """
        Return True if the two objects are equal, False otherwise.

        Objects are considered equal if their `x`, `y`, `z`, `vx`, `vy`, `vz`,
        and `m` attributes are all equal.

        Parameters
        ----------
        other : object
            The other object to compare to.

        Returns
        -------
        bool
        True if the two objects are equal, False otherwise.
        """
        return all(
            [
                np.array_equal(self.x, other.x),
                np.array_equal(self.y, other.y),
                np.array_equal(self.z, other.z),
                np.array_equal(self.vx, other.vx),
                np.array_equal(self.vy, other.vy),
                np.array_equal(self.vz, other.vz),
                np.array_equal(self.m, other.m),
            ]
        )

    def __repr__(self):
        """Representation method.

        Returns
        -------
            str
                Name plus number of points in the box
        """
        cls_name = type(self).__name__
        length = len(self)
        return f"<{cls_name} size={length}>"

    def size(self):
        size = math.ceil(np.max(self.z.value))
        return size
