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
file for VoidFinderTK. It provides functionality to create an empty
configuration, read an existing configuration, and store global settings.

The module automatically creates a default configuration file if it doesn't
exist and loads the settings into a global SETTINGS object.
"""

# =============================================================================
# IMPORTS
# =============================================================================

import datetime as dt
import os
import pathlib

import attrs

import yaml

from . import __version__ as VERSION


# =============================================================================
# CONSTANTS
# =============================================================================

#: pathlib.Path: The default path for the VoidFinderTK configuration file.
DEFAULT_CONF_PATH = pathlib.Path.home() / ".voidfindertk" / "vftk.yaml"

CWD_CONF_PATH = pathlib.Path.cwd() / "vftk.yaml"

# =============================================================================
# CONFIGURATION
# =============================================================================

_RO = {"readonly": True}


@attrs.frozen()
class _Settings:

    _ENV_PREFFIX = "VFTK_"

    voidfindertk_version: str = attrs.field(default=VERSION, metadata=_RO)
    created_at: str = attrs.field(metadata=_RO)
    zobov_path: str = attrs.field(default="")
    popcorn_path: str = attrs.field(default="")

    @created_at.default
    def _created_at_default(self):
        return dt.datetime.now(dt.timezone.utc).isoformat()

    @classmethod
    def from_yaml(cls, buffer):
        data = yaml.safe_load(buffer)
        return cls(**data)

    def update_from_dict(self, data):
        current = self.to_dict(read_only=False)

        diff = set(data).difference(current)
        if diff:
            raise ValueError("Cant assing attribute/s: {diff}")

        current.update(data)
        for k, v in current.items():
            super().__setattr__(k, v)

    def update_from_env(self):
        data = {
            k.replace(self._ENV_PREFFIX, "", 1).lower(): v
            for k, v in os.environ.items()
            if k.startswith(self._ENV_PREFFIX)
        }
        self.update_from_dict(data)

    def update_from_yaml(self, buffer):
        data = self.from_yaml(buffer).to_dict(read_only=False)
        self.update_from_dict(data)

    def to_dict(self, read_only=True):
        def no_privates(a, _):
            return not a.name.startswith("_") and (
                read_only or a.metadata != _RO
            )

        data = attrs.asdict(self, filter=no_privates)
        return data

    def to_yaml(self, buffer=None, **kwargs):
        data = self.to_dict()
        return yaml.safe_dump(data, stream=buffer, **kwargs)


# =============================================================================
# LOAD CONFIGURATION
# =============================================================================

# create the directory if needed
DEFAULT_CONF_PATH.parent.mkdir(parents=True, exist_ok=True)

if DEFAULT_CONF_PATH.exists():
    with open(DEFAULT_CONF_PATH, "r") as fp:
        SETTINGS = _Settings.from_yaml(fp)
else:
    print("Creating new configuration...")
    SETTINGS = _Settings()
    with open(DEFAULT_CONF_PATH, "w") as fp:
        SETTINGS.to_yaml(fp, default_flow_style=False, sort_keys=False)
    print(f"Please configure {DEFAULT_CONF_PATH}")

if CWD_CONF_PATH.exists():
    with open(CWD_CONF_PATH, "r") as fp:
        SETTINGS.update_from_yaml(fp)

SETTINGS.update_from_env()
