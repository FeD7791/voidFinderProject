import pytest

import numpy as np

from voidfindertk.box import Box


@pytest.fixture(scope="session")
def mkbox_params():
    def _maker(
        *,
        seed=None,
        coordinates_scale=500,
        velocity_scale=220,
        mass_scale=12,
        size=1000,
    ):

        rng = np.random.default_rng(seed=seed)
        params = {
            "x": coordinates_scale * rng.random(size=size),
            "y": coordinates_scale * rng.random(size=size),
            "z": coordinates_scale * rng.random(size=size),
            "vx": velocity_scale * rng.random(size=size),
            "vy": velocity_scale * rng.random(size=size),
            "vz": velocity_scale * rng.random(size=size),
            "m": mass_scale * rng.random(size=size),
        }
        return params

    return _maker


@pytest.fixture(scope="session")
def mkbox(mkbox_params):
    def _maker(**kwargs):
        params = mkbox_params(**kwargs)
        return Box(**params)

    return _maker
