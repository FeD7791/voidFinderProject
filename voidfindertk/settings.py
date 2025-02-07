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

"""Configuration management for VoidFinderTK.

This module handles the creation, reading, and management of the configuration
file for VoidFinderTK. It provides functionality to create an empty
configuration, read an existing configuration, and store global settings.

The module automatically creates a default configuration file if it doesn't
exist and loads the settings into a global SETTINGS object.

Notes
-----
The configuration is stored in YAML format and can be modified through:
    - Direct file manipulation
    - Environment variables (prefixed with VFTK_)
    - Local configuration file in the current working directory

"""

# =============================================================================
# IMPORTS
# =============================================================================

import datetime as dt
import logging
import os
import pathlib

import attrs

import yaml

from . import __version__ as VERSION


# =============================================================================
# CONSTANTS
# =============================================================================

#: Path to the default configuration file
DEFAULT_CONF_PATH = pathlib.Path.home() / ".voidfindertk" / "vftk.yaml"

#: Path to the current working directory configuration file
CWD_CONF_PATH = pathlib.Path.cwd() / "vftk.yaml"

#: Readonly metadata for attrs
_RO = {"readonly": True}

# =============================================================================
# SETUP
# =============================================================================

# Create the configuration directory if needed
DEFAULT_CONF_PATH.parent.mkdir(parents=True, exist_ok=True)

# Set up module logger
logger = logging.getLogger(__name__)

# =============================================================================
# SETTINGS
# =============================================================================


@attrs.frozen()
class _Settings:
    """Configuration settings manager for VoidFinderTK.

    This class handles all configuration settings for the VoidFinderTK
    application, including version information, creation timestamp, and paths
    to required external tools.

    Parameters
    ----------
    voidfindertk_version : str
        Current version of VoidFinderTK
    created_at : str
        ISO format timestamp of when the configuration was created
    zobov_path : str
        Path to the ZOBOV executable
    popcorn_path : str
        Path to the POPCORN executable

    Attributes
    ----------
    _ENV_PREFFIX : str
        Prefix used for environment variables to configure settings

    Notes
    -----
    This class is frozen (immutable) and uses attrs for attribute management.
    Configuration can be updated through methods that create new instances.
    """

    _ENV_PREFFIX = "VFTK_"

    voidfindertk_version: str = attrs.field(default=VERSION, metadata=_RO)
    created_at: str = attrs.field(metadata=_RO)
    zobov_path: str = attrs.field(default="")
    popcorn_path: str = attrs.field(default="")

    @created_at.default
    def _created_at_default(self):
        """Generate creation timestamp in UTC."""
        return dt.datetime.now(dt.timezone.utc).isoformat()

    @classmethod
    def from_yaml(cls, buffer):
        """Create settings instance from YAML buffer.

        Parameters
        ----------
        buffer : file-like object
            YAML formatted configuration data

        Returns
        -------
        _Settings
            New settings instance with loaded configuration
        """
        data = yaml.safe_load(buffer)
        return cls(**data)

    def update_from_dict(self, data):
        """Update settings from dictionary.

        Parameters
        ----------
        data : dict
            Dictionary containing new setting values

        Raises
        ------
        ValueError
            If dictionary contains invalid setting names
        """
        current = self.to_dict(read_only=False)

        diff = set(data).difference(current)
        if diff:
            raise ValueError(f"Cannot assign attribute(s): {diff}")

        current.update(data)
        for k, v in current.items():
            super().__setattr__(k, v)

    def update_from_env(self):
        """Update settings from environment variables.

        Looks for variables prefixed with VFTK_ and updates corresponding
        settings.

        """
        data = {
            k.replace(self._ENV_PREFFIX, "", 1).lower(): v
            for k, v in os.environ.items()
            if k.startswith(self._ENV_PREFFIX)
        }
        self.update_from_dict(data)

    def update_from_yaml(self, buffer):
        """Update settings from YAML buffer.

        Parameters
        ----------
        buffer : file-like object
            YAML formatted configuration data
        """
        data = self.from_yaml(buffer).to_dict(read_only=False)
        self.update_from_dict(data)

    def to_dict(self, read_only=True):
        """Convert settings to dictionary.

        Parameters
        ----------
        read_only : bool, optional
            Whether to include read-only attributes, by default True

        Returns
        -------
        dict
            Dictionary representation of settings
        """

        def no_privates(a, _):
            return not a.name.startswith("_") and (
                read_only or a.metadata != _RO
            )

        data = attrs.asdict(self, filter=no_privates)
        return data

    def to_yaml(self, buffer=None, **kwargs):
        """Convert settings to YAML format.

        Parameters
        ----------
        buffer : file-like object, optional
            Buffer to write YAML to, by default None
        **kwargs
            Additional arguments passed to yaml.safe_dump

        Returns
        -------
        str or None
            YAML string if no buffer provided, None otherwise
        """
        data = self.to_dict()
        return yaml.safe_dump(data, stream=buffer, **kwargs)


def load_default_settings():
    """Load default settings.

    Returns
    -------
    _Settings
        Default settings

    Notes
    -----
    If no configuration file is found, it is created with default values.

    If a local configuration file is found, it is loaded and merged with the
    default settings.

    Environment variables are also read and used to update the settings.

    Priority:
        1. Environment variables
        2. Local configuration file
        3. Default configuration

    """
    # Load or create default configuration
    if DEFAULT_CONF_PATH.exists():
        with open(DEFAULT_CONF_PATH, "r") as fp:
            the_settings = _Settings.from_yaml(fp)
    else:
        logger.warning(
            "No configuration file found. Creating new configuration..."
        )
        the_settings = _Settings()
        with open(DEFAULT_CONF_PATH, "w") as fp:
            the_settings.to_yaml(fp, default_flow_style=False, sort_keys=False)
        logger.warning(
            f"Please configure VoidFinderTK at: {DEFAULT_CONF_PATH}"
        )

    # Load local configuration if it exists
    if CWD_CONF_PATH.exists():
        with open(CWD_CONF_PATH, "r") as fp:
            the_settings.update_from_yaml(fp)

    # Update settings from environment variables
    the_settings.update_from_env()

    return the_settings


#: Global settings
SETTINGS = load_default_settings()
