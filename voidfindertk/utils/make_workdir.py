#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023 - 2024, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
#                            Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
"""Module for creation of temporary directories."""
import datetime as dt
import pathlib
import tempfile


def create_run_work_dir(*, workdir_path):
    """
    Create a temporal directory inside the working directory of the ZobovVF\
    class workdir.

    Returns
    -------
        run_work_dir: pathlib.Path
            path of the work directoty
    """
    timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
    run_work_dir = pathlib.Path(
        tempfile.mkdtemp(suffix=timestamp, dir=workdir_path)
    )
    return run_work_dir
