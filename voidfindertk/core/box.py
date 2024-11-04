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
"""Class Box object constructor."""


# =============================================================================
# IMPORTS
# =============================================================================
import math

from astropy import units as u

import attrs

import numpy as np

import uttr

from . import plot_acc


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

    """

    x = uttr.ib(converter=np.array, unit=u.Mpc)
    y = uttr.ib(converter=np.array, unit=u.Mpc)
    z = uttr.ib(converter=np.array, unit=u.Mpc)
    vx = uttr.ib(converter=np.array, unit=u.Mpc / u.second)
    vy = uttr.ib(converter=np.array, unit=u.Mpc / u.second)
    vz = uttr.ib(converter=np.array, unit=u.Mpc / u.second)
    m = uttr.ib(converter=np.array, unit=u.M_sun)

    plot = uttr.ib(
        init=False,
        default=attrs.Factory(plot_acc.BoxPlotter, takes_self=True),
    )

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        box_attributes = [
            self.x.value,
            self.y.value,
            self.z.value,
            self.vx.value,
            self.vy.value,
            self.vz.value,
            self.m.value,
        ]
        # set of number of elements of each array attribute.
        lengths = set(map(len, box_attributes))

        # Validator 1 : if another lenght is found then lengths hast more than
        # one value.
        if len(lengths) != 1:
            raise ValueError("Arrays should be of the same size")

        # Validator 2 : check if the box is a cube.
        box_side = set(
            map(
                lambda arr: math.ceil(np.max(arr)) - math.floor(np.min(arr)),
                box_attributes[:3],
            )
        )

        if len(box_side) != 1:
            raise ValueError(
                "Not a cube: "
                f"xmax: {math.ceil(np.max(self.x.value))} "
                f"ymax: {math.ceil(np.max(self.y.value))} "
                f"zmax: {math.ceil(np.max(self.z.value))} "
                f"xmin: {math.floor(np.min(self.x.value))} "
                f"ymin: {math.floor(np.min(self.y.value))} "
                f"zmin: {math.floor(np.min(self.z.value))} "
            )

    def __len__(self):
        """Length method.

        Returns
        -------
            int
                the number of elements in the box
        """
        return len(self.x)

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
        return f"<{cls_name} size={length} xyzmin={self.min()} \
            xyzmax={self.max()}>"

    def size(self):
        """
        Returns the lenght of side of the box.

        Returns
        -------
            int : Lenght of box.
        """
        return math.ceil(np.max(self.z.value)) - math.floor(
            np.min(self.z.value)
        )

    def min(self):
        """
        Returns the minimun value, in length position of the box.

        This value is the same for the x,y,z coordinates since all tracers are
        aligned in a box.

        Returns
        -------
            int
                The lenght minimun value of the box.
        """
        return math.floor(np.min(self.z.value))

    def max(self):
        """
        Returns the maximun value, in length position of the box.

        This value is the same for the x,y,z coordinates since all tracers are
        aligned in a box.

        Returns
        -------
            int
                The lenght maximun value of the box.
        """
        return math.ceil(np.max(self.z.value))

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
