###############################################################################
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Federico, Gualpa Sebastian, Cabral Juan
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
###############################################################################
import pandas as pd


def box_to_csv(box, file_path, **kwargs):
    kwargs.setdefault("path_or_buf", file_path)
    kwargs.setdefault("delim_whitespace", True)
    kwargs.setdefault("header", False)
    kwargs.setdefault("index", False)
    df = pd.DataFrame(box.__dict__)
    df.to_csv(**kwargs)
