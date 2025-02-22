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


def _box_converter(a):
    """
    Logic for attrs/uttrs converter.

    Parameters
    ----------
        a : array or None
    Returns
    -------
        a : np.array
    """
    if a is not None:
        a = np.array(a)
    else:
        a = np.array([])
    return a


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

    x = uttr.ib(default=None, converter=_box_converter, unit=u.Mpc)
    y = uttr.ib(default=None, converter=_box_converter, unit=u.Mpc)
    z = uttr.ib(default=None, converter=_box_converter, unit=u.Mpc)
    vx = uttr.ib(default=None, converter=_box_converter, unit=u.Mpc / u.second)
    vy = uttr.ib(default=None, converter=_box_converter, unit=u.Mpc / u.second)
    vz = uttr.ib(default=None, converter=_box_converter, unit=u.Mpc / u.second)
    m = uttr.ib(
        default=None,
        converter=_box_converter,
    )

    plot = uttr.ib(
        init=False,
        default=attrs.Factory(plot_acc.BoxPlotter, takes_self=True),
    )

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        attributes = np.array([self.x, self.y, self.z], dtype=object)
        k = np.array(list(map(len, attributes)))
        indx = np.where(k > 0)[0]
        lengths = set(map(len, attributes[indx]))

        # Validator 1 : if another lenght is found then lengths hast more than
        # one value.
        if len(lengths) != 1:
            raise ValueError("Arrays should be of the same size")

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
        return f"<{cls_name} size={length}"

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
        return int(
            np.min(
                np.array(
                    [
                        math.floor(np.min(self.x.value)),
                        math.floor(np.min(self.y.value)),
                        math.floor(np.min(self.z.value)),
                    ]
                )
            )
        )

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
        return int(
            np.max(
                np.array(
                    [
                        math.ceil(np.max(self.x.value)),
                        math.ceil(np.max(self.y.value)),
                        math.ceil(np.max(self.z.value)),
                    ]
                )
            )
        )

    def mass_cutoff(self, mass_threshold):
        """
        Filter the points by mass threshold.

        Returns a new `Box` object that only includes points with a mass
        greater than the given threshold.

        Parameters
        ----------
        mass_threshold : float
            The mass threshold above which points will be included.

        Returns
        -------
        Box
            A new `Box` object containing only points with mass greater
            than the threshold.
        """
        idx = np.where(self.m > mass_threshold)[0]
        properties = [
            self.x,
            self.y,
            self.z,
            self.vx,
            self.vy,
            self.vz,
            self.m,
        ]
        names = ["x", "y", "z", "vx", "vy", "vz", "m"]
        available_properties = np.array(list(map(len, properties)))

        new_properties = {}
        for e in zip(names, available_properties != 0, properties):
            if e[1]:
                new_properties[e[0]] = e[2][idx]
            else:
                new_properties[e[0]] = np.array([])
        return self.copy(**new_properties)

    def to_dict(self):
        """Method used to convert the class Box to a dictionary.

        Returns
        -------
            dict
        """
        return attrs.asdict(self, filter=lambda a, _: a.init)

    def copy(self, **kwargs):
        """Method used to perform a deep copy of the class Box.

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments to update the Box parameters.

        Retruns
        -------
            Box object with a copy of Box Parameters.
        """
        cls = type(self)

        data = self.to_dict()
        data.update(kwargs)

        return cls(**data)
