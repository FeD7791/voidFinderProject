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

"""Module for process data obtained from PopCorn void finder."""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np


# Popcorn reader (Old version - to be changed)
def read_pop(filename):
    """
    Gets output from PopCorn Voidfinder from its output file.

    Reads a data table of popcorn voids from an ASCII file (without header)
    and structures it as a dictionary.

    Parameters
    ----------
    filename : str
        The path to the data file.

    Returns
    -------
    popcorn : dict
        A dictionary with keys corresponding to popcorn voids attributes.

    Popcorn voids attributes
    -------------------------
    id : int
        Popcorn ID.
    nmem : int
        Number of member spheres.
    vol : float
        Popcorn volume.
    reff : float
        Effective radius.
    npart : int
        Number of particles inside.
    flag : int
        For internal control, can be ignored.
    pop : dict
        A dictionary with attributes of the member spheres (length: nmem).
        * x : list of float
            x-coordinates of the centres of the member spheres.
        * y : list of float
            y-coordinates of the member spheres.
        * z : list of float
            z-coordinates of the member spheres.
        * r : list of float
            Radii of the member spheres.
        * fvol : list of float
            Volume contribution to the total volume.
        * level : list of int
            Hierarchy level: 0 for main sphere, 1 for secondary spheres, etc.
    """
    with open(filename, "r") as file:
        npop = int(
            file.readline().strip()
        )  # number of popcorn voids (integer)
        print(f"Number of popcorns: {npop}")

        popcorn = {
            "id": [],
            "nmem": [],
            "vol": [],
            "reff": [],
            "npart": [],
            "flag": [],
            "pop": [],
        }

        for p in range(npop):
            head_pop = list(
                file.readline().strip().split()
            )  # this line is the header
            id = int(head_pop[0])
            nmem = int(head_pop[1])
            vol = float(head_pop[2])
            npart = int(head_pop[3])
            flag = int(head_pop[4])

            reff = (vol * 3 / (4 * np.pi)) ** (1.0 / 3.0)
            # NOTE: this is a derived result, not contained in the original
            # catalogue

            # Building the `pop` object:
            x = []
            y = []
            z = []
            r = []
            fvol = []
            level = []
            if nmem > 0:
                for _ in range(nmem):
                    popmem = list(
                        file.readline().strip().split()
                    )  # attributes of the member spheres
                    x.append(float(popmem[0]))
                    y.append(float(popmem[1]))
                    z.append(float(popmem[2]))
                    r.append(float(popmem[3]))
                    fvol.append(float(popmem[4]))
                    level.append(int(popmem[5]))
            pop = {
                "x": x,
                "y": y,
                "z": z,
                "r": r,
                "fvol": fvol,
                "level": level,
            }  # this is an individual popcorn

            if npart > 0:
                for _ in range(npart):
                    file.readline()  # reading inner particles ID (ignoring)

            popcorn["id"].append(id)
            popcorn["nmem"].append(nmem)
            popcorn["vol"].append(vol)
            popcorn["reff"].append(reff)
            popcorn["npart"].append(npart)
            popcorn["flag"].append(flag)
            popcorn["pop"].append(pop)

    return popcorn