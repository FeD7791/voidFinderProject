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
"""Module for interacting with the SVF Popcorn void finder and related\
utilities."""
# =============================================================================
# IMPORTS
# =============================================================================
from configparser import ConfigParser

import pandas as pd

import sh

from ..utils import chdir


def config_file_maker(
    *,
    # Parameters
    trsfile,
    filefmt,
    num_file,
    sphfile,
    popfile,
    auxfiles,
    rawpopfile,
    pairsfile,
    boxsize,
    densth,
    minradius,
    maxradius,
    massmin,
    eps,
    # Directory path to place the file
    path,
):
    """
    Builds configuration file with void run parameters.

    Generates the configuration file with the parameters to run the void
    finder in the desired path.

    Parameters
    ----------
    trsfile : str
        Input tracer file.
    filefmt : str
        File format options: "ASCII", "STREAM", "HDF5",
        "HDF5_SUBFIND_GROUPS", "HDF5_SUBFIND_SUBHALOS", "GADGET1",
        "GADGET2", "GADGET4_TYPE1".
    num_file : str
        Number of files for the tracer catalogue, especially important for
        Gadget outputs.
    sphfile : str
        Output Spherical voids catalogue.
    popfile : str
        Output Popcorn voids catalogue (after cleaning overlaps).
    auxfiles : str
        Whether to use auxiliary files: {"true", "false"}.
    rawpopfile : str
        File with Popcorn voids before cleaning overlaps.
    pairsfile : str
        File with pairs of touching Popcorn voids.
    boxsize : str
        Length of the box in the same units as the tracer input file
        (required only for ASCII inputs, otherwise ignored).
    densth : str
        Integrated density threshold for identification.
    minradius : str
        Minimal radii allowed for a sphere member in input units.
    maxradius : str
        Maximal radii allowed for a sphere member in input units.
    massmin : str
        Minimal halo mass allowed (set to 0 if not applicable).
    eps : str
        Obsolete flag; do not modify.
    path : pathlib.Path
        Path where the file will be generated.
    """
    config = ConfigParser(allow_no_value=True)
    config.optionxform = str
    config.add_section("INPUT_PARAMS")
    config.set("INPUT_PARAMS", "TRSFILE", trsfile)
    config.set("INPUT_PARAMS", "FILEFMT", filefmt)
    config.set("INPUT_PARAMS", "NUM_FILE", num_file)
    config.set("INPUT_PARAMS", "SPHFILE", sphfile)
    config.set("INPUT_PARAMS", "POPFILE", popfile)
    config.set("INPUT_PARAMS", "AUXFILES", auxfiles)
    config.set("INPUT_PARAMS", "RAWPOPFILE", rawpopfile)
    config.set("INPUT_PARAMS", "PAIRSFILE", pairsfile)
    config.set("INPUT_PARAMS", "BOXSIZE", boxsize)
    config.set("INPUT_PARAMS", "DENSTH", densth)
    config.set("INPUT_PARAMS", "MINRADIUS", minradius)
    config.set("INPUT_PARAMS", "MAXRADIUS", maxradius)
    config.set("INPUT_PARAMS", "MASSMIN", massmin)
    config.set("INPUT_PARAMS", "EPS", eps)

    with open(path, "w") as configfile:
        config.write(configfile)


def popcorn_svf_input_data_builder(*, box, file_path):
    """
    Generates input file from box in the desired file path.

    Parameters
    ----------
    box : Box
        Box object containing the tracers and their properties.
    file_path : pathlib.Path
        File path to place the generated input file.
    """
    df = pd.DataFrame(box.__dict__)
    # df.drop(labels=["_len"], axis=1, inplace=True)
    # Popcorn input file format
    df = df[["m", "x", "y", "z", "vx", "vy", "vz"]]
    df.to_csv(
        file_path, sep=" ", index=False, header=False, float_format="%.2f"
    )


def spherical_popcorn_void_finder(
    *, mpi_flags, bin_path, conf_file_path, work_dir_path
):
    """
    Runs the Popcorn-Spherical Void Finder using MPI.

    Parameters
    ----------
    mpi_flags : str
        MPI flags for execution.
    bin_path : pathlib.Path
        Path to the folder containing the binary file for the void finder.
    conf_file_path : pathlib.Path
        Path to the configuration file used to run the void finder.
    work_dir_path : str
        Path to the working directory.

    Returns
    -------
    output : str
        Standard output from the run.
    """
    # Reference to mpi
    # https://docs.oracle.com/cd/E19356-01/820-3176-10/ExecutingPrograms.html
    svf_mpi = sh.Command("svf", search_paths=[bin_path])
    params = "config=" + str(conf_file_path)
    # Command will be executed from work_dir_path path.
    with chdir(work_dir_path):
        output = svf_mpi(params)
    return output
    # subprocess.run(["mpirun", "-np", "1","--bind-to","none",
    # "./svf", "./configuration/vars.conf"])
