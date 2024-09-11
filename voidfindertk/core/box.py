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


class DataBox:
    box = None

    def __init__(self, aux_box):
        self.box = aux_box


@uttr.s(repr=False, frozen=True, cmp=False)
class Box:
    """Box Class.

    Class used to describe a set of points (x,y,z) alongside with its
    velocities (vx,vy,vz)

    Attributes
    ----------
    x : numpy.ndarray
        (x,y,z) array of position elements
    y : numpy.ndarray
        (x,y,z) array of position elements
    z : numpy.ndarray
        (x,y,z) array of position elements
    vx : numpy.ndarray
        (vx,vy,vz) array of velocity elements
    vy : numpy.ndarray
        (vx,vy,vz) array of velocity elements
    vz : numpy.ndarray
        (vx,vy,vz) array of velocity elements
    m : numpy.ndarray
        array of masses

    Methods
    -------
    build_df()
        Transforms the set into a pandas dataframe

    _check_length()
        Verifies that the length of the inputs are the same

    slice(n, parameter)
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

        for e in (self.x, self.y, self.z, self.vx, self.vy, self.vz, self.m):
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
                "Not a cube: "
                f"xmax: {math.ceil(np.max(self.x.value))} "
                f"ymax: {math.ceil(np.max(self.y.value))} "
                f"zmax: {math.ceil(np.max(self.z.value))}"
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
        """
        Representation method.

        Returns
        -------
            str
                Name plus number of points in the box
        """
        cls_name = type(self).__name__
        length = len(self)
        return f"<{cls_name} size={length}>"

    def size(self):
        """
        Method used to determine the size of the box along the z-axis.

        Returns
        -------
        int
            The size of the box along a side box-axis
        """
        size = math.ceil(np.max(self.z.value))
        return size

    def copy(self):
        """
        Method used to perform a deep copy of the class Box.

        Retruns
        -------
            Box object with a copy of Box Parameters.
        """
        cls = type(self)
        new_tracers = cls(
            x=self.x,
            y=self.y,
            z=self.z,
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            m=self.m,
        )
        return new_tracers
