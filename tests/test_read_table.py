import numpy as np

import pandas as pd

import pytest

from voidfindertk import box, read_table


def make_dataset(xyz_size=500, vxyz_size=220, m_size=12, n=1000, random_seed=42):
    rng = np.random.default_rng(seed=random_seed)
    dataset = pd.DataFrame(
        {
            "x": xyz_size * rng.random(n),
            "y": xyz_size * rng.random(n),
            "z": xyz_size * rng.random(n),
            "vx": vxyz_size * rng.random(n),
            "vy": vxyz_size * rng.random(n),
            "vz": vxyz_size * rng.random(n),
            "m": m_size * rng.random(n),
        }
    )
    with open("./tests/datasets/data.txt", "w") as f:
        dataset.to_string(f, col_space=10)


def test_read_table_type_output():
    assert isinstance(read_table("./tests/datasets/data_random_1.txt"), box.Box)


def test_nullvalues_on_read():
    with pytest.raises(TypeError):
        read_table("./tests/datasets/data_random_2.txt")


def test_number_of_columns():
    with pytest.raises(ValueError):
        read_table("./tests/datasets/data_random_1.txt", usecols=[1, 2, 3])
