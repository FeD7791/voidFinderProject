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

"""Module for interacting with the Popcorn void finder and related \
utilities."""

# =============================================================================
# IMPORTS
# =============================================================================

import configparser

import sh

from ..utils import chdir

# =============================================================================
# FUNCTIONS
# =============================================================================

# Reference:
# https://gitlab.com/dante.paz/popcorn_void_finder#43-popcorn-void-finder


def popcorn_void_finder(*, mpi_flags, bin_path, conf_file_path, work_dir_path):
    """
    Executes the Popcorn void finder with the specified configuration.

    Parameters
    ----------
    mpi_flags : str
        Flags for MPI configuration (not currently used in the command).
    bin_path : str
        Path to the directory containing the Popcorn executable.
    conf_file_path : str
        Path to the configuration file to be used by the Popcorn void finder.
    work_dir_path : str
        Directory path where the command will be executed.

    Returns
    -------
    output : str
        The output from the Popcorn command execution.
    """
    popcorn = sh.Command("popcorn", search_paths=[bin_path])
    params = "config=" + str(conf_file_path)
    # Command will be executed from work_dir_path path.
    with chdir(work_dir_path):
        output = popcorn(params)
    return output


def compute_intersects(*, bin_path, conf_file_path, work_dir_path):
    """
    Executes the compute_intersecs command with the specified configuration.

    Parameters
    ----------
    bin_path : str
        Path to the directory containing the compute_intersecs executable.
    conf_file_path : str
        Path to the configuration file to be used by the compute_intersecs
        command.
    work_dir_path : str
        Directory path where the command will be executed.

    Returns
    -------
    output : str
        The output from the compute_intersecs command execution.
    """
    compute_intersecs = sh.Command(
        "compute_intersecs", search_paths=[bin_path]
    )
    params = "config=" + str(conf_file_path)
    # Command will be executed from work_dir_path path.
    with chdir(work_dir_path):
        output = compute_intersecs(params)
    return output


def clean_duplicates(*, bin_path, conf_file_path, work_dir_path):
    """
    Executes the clean_duplicates command with the specified configuration.

    Parameters
    ----------
    bin_path : str
        Path to the directory containing the clean_duplicates executable.
    conf_file_path : str
        Path to the configuration file to be used by the clean_duplicates
        command.
    work_dir_path : str
        Directory path where the command will be executed.

    Returns
    -------
    output : str
        The output from the clean_duplicates command execution.
    """
    clean_duplicates = sh.Command("clean_duplicates", search_paths=[bin_path])
    params = "config=" + str(conf_file_path)
    # Command will be executed from work_dir_path path.
    with chdir(work_dir_path):
        output = clean_duplicates(params)
    return output


def read_and_modify_config(*, config_file_path, section, parameter, new_value):
    """
    Modifies a specified parameter in a configuration file.

    Parameters
    ----------
    config_file_path : str
        Path to the configuration file to be modified.
    section : str
        The section in the configuration file that contains the parameter to
        modify.
    parameter : str
        The parameter within the section that needs to be updated.
    new_value : str
        The new value to set for the specified parameter.

    Returns
    -------
    None
        This function modifies the configuration file in place and does not
        return a value.
    """
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    config.optionxform = str
    # Read the configuration file
    config.read(config_file_path)

    # Check if the section exists
    if not config.has_section(section):
        print(f"Section '{section}' not found in the configuration file.")
        return

    # Check if the parameter exists
    if not config.has_option(section, parameter):
        print(f"Parameter '{parameter}' not found in section '{section}'.")
        return

    # Modify the parameter value
    config.set(section, parameter, new_value)

    # Save the changes back to the configuration file
    with open(config_file_path, "w") as configfile:
        config.write(configfile)
