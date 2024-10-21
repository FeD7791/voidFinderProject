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

"""
Configuration management for VoidFinderTK.

This module handles the creation, reading, and management of the configuration
file for VoidFinderTK. It provides functionality to create an empty configuration,
read an existing configuration, and store global settings.

The module automatically creates a default configuration file if it doesn't exist
and loads the settings into a global SETTINGS object.
"""

# =============================================================================
# IMPORTS
# =============================================================================

import datetime as dt
import json
import os
import pathlib

from . import __version__ as VERSION
from .utils.bunch import Bunch

# =============================================================================
# DEFAULT CONF
# =============================================================================

#: dict: The default empty configuration structure.
_EMPTY_CONF = {
    "void_finder_tk_version": VERSION,
    "created_at": None,
    "paths": {
        "zobov_path": "",
        "popcorn_path": "",
    },
}

# =============================================================================
# API
# =============================================================================


def create_empty_conf(fp):
    """
    Create an empty configuration and write it to the given file object.

    Parameters
    ----------
    fp : file object
        The file object to write the configuration to.

    Notes
    -----
    The created configuration includes the current version of VoidFinderTK
    and the creation timestamp in UTC.
    """
    conf = _EMPTY_CONF.copy()
    conf["created_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    json.dump(conf, fp, indent=2)


def read_conf(fp):
    """
    Read a configuration from the given file object.

    Parameters
    ----------
    fp : file object
        The file object to read the configuration from.

    Returns
    -------
    Bunch
        A Bunch object containing the configuration data.
    """
    return Bunch(fp.name, json.load(fp))


# =============================================================================
# CONSTANTS
# =============================================================================

#: pathlib.Path: The path to the user's home directory.
USER_HOME_PATH = pathlib.Path(os.path.expanduser("~"))

#: pathlib.Path: The default path for the VoidFinderTK configuration file.
DEFAULT_CONF_PATH = USER_HOME_PATH / ".voidfindertk.json"

# Create a new configuration file if it doesn't exist
if not DEFAULT_CONF_PATH.exists():
    print("Creating new configuration...")
    with open(DEFAULT_CONF_PATH, "w") as fp:
        create_empty_conf(fp)
    print(f"Please configure {DEFAULT_CONF_PATH}")

# Load the settings from the configuration file
with open(DEFAULT_CONF_PATH) as fp:
    #: Bunch: Global settings object loaded from the configuration file.
    SETTINGS = read_conf(fp)
