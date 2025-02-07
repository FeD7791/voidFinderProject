#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2023 - 2024, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
#                            Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.

# =============================================================================
# IMPORTS
# =============================================================================

import os
import pathlib
import tempfile

import numpy as np

from voidfindertk.popcorn import _pc_postprocessing

# =============================================================================
# TESTS
# =============================================================================


def test_get_properties():
    data = """40210
    108 1 1.987356e+03 1 -1
    3.446931e+02 4.357584e+02 1.575218e-01 7.799421e+00 1.987356e+03 0
    23378
    129 1 3.008339e+03 1 -1
    8.203903e+02 5.071564e+02 3.739926e+00 8.955284e+00 3.008339e+03 0
    80740
    131 1 3.449483e+03 1 -1
    4.717176e+02 5.122758e+02 3.188705e+00 9.373213e+00 3.449483e+03 0
    27728
    132 1 5.268928e+03 2 -1
    8.431208e+02 5.190329e+02 1.967021e+00 1.079472e+01 5.268928e+03 0
    28178
    28179
    133 1 6.786148e+03 3 -1
    8.886003e+02 5.165095e+02 2.250668e+00 1.174478e+01 6.786148e+03 0
    80799
    6761039
    6761451
    44997 2 9.015377e+04 61 -1
    7.467082e+01 9.542412e+02 7.565505e+02 2.435306e+01 6.049914e+04 0
    8.092774e+01 9.775796e+02 7.595909e+02 2.195427e+01 9.015377e+04 1
    4971569
    5024910
    5024911
    5025418
    5078175
    5078176
    5078632
    5078633
    5078634
    5078635
    5079067
    5079068
    5079069
    5079070
    5079884
    5080276
    5080645
    5081014
    5132058
    5132964
    5135351
    5135352
    5135353
    5187087
    5187488
    5188013
    5188014
    5188414
    5188415
    5188854
    5190174
    5240593
    5240594
    5240595
    5241060
    5241061
    5241062
    5241063
    5241064
    5241065
    5241066
    5241526
    5241527
    5241528
    5242007
    5242410
    5242411
    5242412
    5242810
    5242811
    5242812
    5243204
    5294071
    5294072
    5294521
    5294526
    5295002
    5295836
    5295837
    5295838
    5295839"""

    with tempfile.NamedTemporaryFile(delete=False, mode="r+") as temp_file:
        # Write the data
        temp_file.write(data)

        # path to file
        temp_filename = pathlib.Path(os.path.abspath(temp_file.name))

        # Read the file
        # voids, sphere = _pc_postprocessing.get_properties(
        #     filename=temp_filename
        #     )

        # The format is the following:

        # number of Popvoids
        # ID | n_mem | vol | n_part | flag
        # x | y | z | r | fracvol | level
        # ID-1 (id particle 1 in void)
        # ID-2
        # ...
        # ID-n

    voids, sphere = _pc_postprocessing.get_properties(filename=temp_filename)
    os.remove(temp_filename)

    ids = np.array([108, 129, 131, 132, 133, 44997])

    n_mem = np.array([1, 1, 1, 1, 1, 2])

    vol = np.array(
        [
            1.987356e03,
            3.008339e03,
            3.449483e03,
            5.268928e03,
            6.786148e03,
            9.015377e04,
        ]
    )

    n_part = np.array([1, 1, 1, 2, 3, 61])

    x = np.array(
        [
            3.446931e02,
            8.203903e02,
            4.717176e02,
            8.431208e02,
            8.886003e02,
            7.467082e01,
            8.092774e01,
        ]
    )
    y = np.array(
        [
            4.357584e02,
            5.071564e02,
            5.122758e02,
            5.190329e02,
            5.165095e02,
            9.542412e02,
            9.775796e02,
        ]
    )
    z = np.array(
        [
            1.575218e-01,
            3.739926e00,
            3.188705e00,
            1.967021e00,
            2.250668e00,
            7.565505e02,
            7.595909e02,
        ]
    )

    r = np.array(
        [
            7.799421e00,
            8.955284e00,
            9.373213e00,
            1.079472e01,
            1.174478e01,
            2.435306e01,
            2.195427e01,
        ]
    )

    level = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
    assert np.allclose(np.array(sphere.x, dtype=float), x, atol=0.001)
    assert np.allclose(np.array(sphere.y, dtype=float), y, atol=0.001)
    assert np.allclose(np.array(sphere.z, dtype=float), z, atol=0.001)
    assert np.allclose(np.array(sphere.radius, dtype=float), r, atol=0.001)

    assert np.allclose(np.array(voids.id, dtype=int), ids, atol=0.1)
    assert np.allclose(np.array(voids.n_mem, dtype=int), n_mem, atol=0.1)
    assert np.allclose(np.array(voids.n_part, dtype=int), n_part, atol=0.1)

    assert np.allclose(np.array(voids.volume, dtype=float), vol, atol=0.001)
    assert np.allclose(np.array(sphere.level, dtype=float), level, atol=0.1)
