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

"""Global configuration for pytest suite.

"""

# =============================================================================
# IMPORTS
# =============================================================================

import io
import os
import pathlib

import grispy as gsp

import joblib

import numpy as np

import pandas as pd

import pytest

from voidfindertk.core.box import Box

# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))

MOCK_DATA = PATH / "mock_data"

# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture(scope="session")
def load_mock_data():
    def loader(*folder_and_filename):
        full_path = MOCK_DATA.joinpath(*folder_and_filename)
        return joblib.load(full_path)

    return loader


@pytest.fixture(scope="session")
def mkbox_params():  # aca vas a retornar la funcion maker que la podes
    # modificar
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
            "x": rng.uniform(0, coordinates_scale, size=size),
            "y": rng.uniform(0, coordinates_scale, size=size),
            "z": rng.uniform(0, coordinates_scale, size=size),
            "vx": rng.uniform(0, velocity_scale, size=size),
            "vy": rng.uniform(0, velocity_scale, size=size),
            "vz": rng.uniform(0, velocity_scale, size=size),
            "m": rng.uniform(0, mass_scale, size=size),
        }
        return params

    return _maker


@pytest.fixture(scope="session")
def mkbox(mkbox_params):
    def _maker(**kwargs):
        params = mkbox_params(**kwargs)
        return Box(**params)

    return _maker


@pytest.fixture(scope="session")
def random_buffer():
    def _maker(*, empty_row=False, sp_characters=None):
        data = np.random.random((100, 7))
        df = pd.DataFrame(data)
        n = 0
        if empty_row:
            df.loc[100] = [np.nan] * 7
            n = 1
        elif sp_characters:
            for i in np.arange(len(sp_characters)):
                df.loc[100 + n + i] = [sp_characters[i]] * 7
        src = df.to_csv(
            sep=" ",
            index=False,
            float_format="%.5f",
            header=False,
            na_rep=np.nan,
        )
        buff = io.StringIO(src)
        return buff

    return _maker


@pytest.fixture
def make_spherical_voids_params():
    # aca vas a retornar un valor, no una funcion
    def _maker(**kwargs):
        params = {
            "n_voids": 1000,
            "rad_scale": 15,
            "xyz_void_max_scale": 500,
            "vel_xyz_void_max_scale": 200,
            "min_delta": -0.95,
            "max_delta": -0.90,
            "min_poisson": -0.5,
            "max_poisson": 0.5,
            "dtype": 0.5,
            "nran": 200,
            "seed": 42,
        }
        for key, value in kwargs.items():
            params[key] = value

        rng = np.random.default_rng(seed=params["seed"])
        void_params = {
            "rad": rng.uniform(0, params["rad_scale"], params["n_voids"]),
            "x_void": rng.uniform(
                0, params["xyz_void_max_scale"], params["n_voids"]
            ),
            "y_void": rng.uniform(
                0, params["xyz_void_max_scale"], params["n_voids"]
            ),
            "z_void": rng.uniform(
                0, params["xyz_void_max_scale"], params["n_voids"]
            ),
            "vel_x_void": rng.uniform(
                0, params["vel_xyz_void_max_scale"], params["n_voids"]
            ),
            "vel_y_void": rng.uniform(
                0, params["vel_xyz_void_max_scale"], params["n_voids"]
            ),
            "vel_z_void": rng.uniform(
                0, params["vel_xyz_void_max_scale"], params["n_voids"]
            ),
            "delta": rng.uniform(
                params["min_delta"], params["max_delta"], params["n_voids"]
            ),
            "poisson": rng.uniform(
                params["min_poisson"], params["max_poisson"], params["n_voids"]
            ),
            "dtype": rng.uniform(0, params["dtype"], params["n_voids"]),
            "nran": rng.uniform(0, params["nran"], params["n_voids"]),
        }
        return void_params

    return _maker


@pytest.fixture
def build_cloud():
    def _build_cloud(*, seed=2, lmin=0, lmax=1000, n_points=100**3):
        """
        Builds a cloud of n_points , where: lmin < x,y,z < lmax.

        Parameters
        ----------
            seed: int
                Seed for random number generator.
            lmin: float
                Low limit of random generated points.
            lmax: float
                Max limit of random generated points.
            n_points: int
                Number of points to be generated
        Returns
        -------
            cloud: array
                Array of n_points with x,y,z coordinates between lmin and lmax.

        """
        np.random.seed(seed)
        # Create point cloud
        cloud = np.random.uniform(lmin, lmax, size=(n_points, 3))
        return cloud

    return _build_cloud


@pytest.fixture
def build_spherical_void():
    def _build_spherical_void(
        delta,
        centers: np.ndarray,
        radii: float,
        cloud,
    ):
        """
        Takes a cloud of tracers and removes tracers inside an spherical shell
        around each void center so that each void has a desired density
        density contrast.

        Parameters
        ----------
            delta: float
                Integrated density contrast of the void. -1 < delta < 0
            centers: array
                Array of x,y,z positions of void centers.
            radii: float
                Radius of void. (All void have same radius)
            cloud: array
                Collection of tracers, each with positions x,y,z that constitu-
                te the universe (Box)
        Returns
        -------
            cloud_with_voids : array
                Collection of tracers of cloud input, minus some tracers around
                each void so that these have the desired density contrast.
        """
        cloud_volume = round(np.max(cloud) - np.min(cloud)) ** 3
        cloud_density = len(cloud) / cloud_volume
        density_voids = (1 + delta) * cloud_density

        # Build Grid
        grid = gsp.GriSPy(cloud)
        # Set periodicity
        periodic = {
            0: (np.round(np.min(cloud), 0), np.round(np.max(cloud), 0)),
            1: (np.round(np.min(cloud), 0), np.round(np.max(cloud), 0)),
            2: (np.round(np.min(cloud), 0), np.round(np.max(cloud), 0)),
        }
        grid.set_periodicity(periodic, inplace=True)

        # Find nearest neighbors
        dist, index = grid.bubble_neighbors(
            centers, distance_upper_bound=radii
        )

        # Calculate right number of tracers so the void gets at the desired
        # density
        n = int(round(density_voids * (4 / 3) * np.pi * (radii**3), 0))

        new_index = [i[: len(i) - int(n - 1)] for i in index]
        # dens_voids

        # Remove from cloud the indicated indexes
        mask = np.ones(len(cloud), dtype=bool)
        for i in new_index:
            mask[i] = False

        cloud_with_voids = cloud[mask]
        # Check the voids have the correct density
        dens_validator(
            cloud_with_voids=cloud_with_voids,
            centers_of_voids=centers,
            radius=radii,
            density=density_voids,
        )
        return cloud_with_voids

    return _build_spherical_void


# =============================================================================
# OTHER
# =============================================================================
def dens_validator(*, cloud_with_voids, centers_of_voids, radius, density):
    """Validator for build_spherical_void. Checks if the created cloud with
    voids has the desired density contrast for each void.
    """
    # Build Grid
    grid = gsp.GriSPy(cloud_with_voids)
    # Add periodicity
    periodic = {
        0: (np.min(cloud_with_voids), np.max(cloud_with_voids)),
        1: (np.min(cloud_with_voids), np.max(cloud_with_voids)),
        2: (np.min(cloud_with_voids), np.max(cloud_with_voids)),
    }
    grid.set_periodicity(periodic, inplace=True)
    # Perform buble search of neighbors around center
    distances, tracers = grid.bubble_neighbors(centers_of_voids, radius)
    # Calculate densities
    densities = np.array(list(map(len, tracers))) / (
        (4 / 3) * np.pi * radius**3
    )
    # Check densities are the desired density
    check_dens = np.all(
        np.less_equal(densities, density * np.ones(densities.shape))
    )
    if not check_dens:
        raise ValueError("Wrong density of voids")
