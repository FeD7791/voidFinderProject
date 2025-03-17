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

"""Global configuration for pytest suite."""

# =============================================================================
# IMPORTS
# =============================================================================

import io
import os
import pathlib
from functools import reduce


import attr

import grispy as gsp

import joblib

import numpy as np

import pandas as pd

import pytest

from scipy.spatial.distance import cdist


from voidfindertk.core.box import Box
from voidfindertk.core.voids import Voids
from voidfindertk.datasets import spherical_cloud
from voidfindertk.settings import SETTINGS
from voidfindertk.utils import box_to_grid
from voidfindertk.zobov import ZobovVF
from voidfindertk.zobov._zobov import _Paths

# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))

MOCK_DATA = PATH / "mock_data"

# =============================================================================
# OTHER
# =============================================================================


@attr.define
class ZobovElements:
    """Contains Paths, Names and executable names of ZOBOV"""

    OUTPUT_VOZINIT = "output_vozinit"
    OUTPUT_JOZOV_VOIDS = "output_txt"
    PARTICLES_IN_ZONES = "part_vs_zone"
    ZONES_IN_VOID = "zones_vs_voids"
    TRACERS_RAW = "tracers_zobov.raw"
    TRACERS_TXT = "tracers_zobov.txt"
    PARTICLES_VS_ZONES_RAW = f"{PARTICLES_IN_ZONES}.dat"
    PARTICLES_VS_ZONES_ASCII = f"{PARTICLES_IN_ZONES}_ascii.txt"
    OUTPUT_JOZOV_VOIDS_DAT = f"{OUTPUT_JOZOV_VOIDS}.dat"
    ZONES_VS_VOID_RAW = f"{ZONES_IN_VOID}.dat"
    ZONES_VS_VOID_ASCII = f"{ZONES_IN_VOID}_ascii.txt"
    ZOBOV_LOADER_BIN = "zobov_loader.so"
    TRACERS_IN_ZONES_BIN = "tracers_in_zones.so"
    ZONES_IN_VOIDS_BIN = "zones_in_void.so"
    ZOBOV = pathlib.Path(SETTINGS.zobov_path)
    CURRENT_FILE_PATH = _Paths.CURRENT_FILE_PATH


@attr.define
class SVFPopcornElements:
    SVF = pathlib.Path(SETTINGS.popcorn_path)
    CONFFILE = SVF / "configuration"
    CONFIG = "vars.conf"
    TRSFILE = "trsfile.dat"
    SPHFILE = "sphfile.dat"
    POPFILE = "popfile.dat"
    RAWPOPFILE = "rawpopfile.dat"
    PAIRSFILE = "pairsfile.dat"


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
def build_box_with_eq_voids():
    """This will return a box of tracers that include void of some defined
    radii and that are equidistant between them.
    """

    def _maker(*, rad=30, cloud=None, delta=-0.9, sep=10, **mass_kwargs):
        """
        Parameters
        ----------
            rad : float
                Radius of the voids.
            cloud : np.array
                Array of tracers.
            delta = : float
                Integrated density contrast.
        Return
        ------
            b : object
                Box with the cloud with voids.
            threshold : float
                Threshold of the density, that divides the underdensities from
                the overdensities.
            centers : array
                Array with the centers of the voids.
            cloud : array
                Updated cloud with tracers.
        """
        lmax = np.max(cloud)
        # Generate centers
        xyz_ = np.arange(0.0, lmax, rad + sep)
        centers = np.column_stack((xyz_, xyz_, xyz_))
        # Get cloud with spherical voids
        cloud_with_voids = spherical_cloud.build_spherical_void(
            delta=delta, centers=centers, radii=rad, cloud=cloud
        )
        # Get cloud with mass
        cloud_with_voids_mass = spherical_cloud.add_mass_to_cloud(
            cloud=cloud_with_voids, **mass_kwargs
        )
        # Density threshold
        threshold = (
            (1 + delta)
            * len(cloud_with_voids)
            / (np.max(cloud_with_voids) ** 3)
        )
        # Transform the cloud to a box
        # box 1
        cloud_ = cloud_with_voids_mass.T
        x = cloud_[0]
        y = cloud_[1]
        z = cloud_[2]
        m = cloud_[3]
        b = Box(x=x, y=y, z=z, vx=x, vy=y, vz=z, m=m)
        return b, threshold, centers, cloud_with_voids

    return _maker


@pytest.fixture
def build_box_from_cloud():
    def _maker(*, cloud=None):
        x, y, z = np.hsplit(cloud, 3)
        b = Box(
            x=x,
            y=y,
            z=z,
            # I don't care about this parameters.
            # But i have to fill them.
            vx=x,
            vy=y,
            vz=z,
            m=z,
        )
        return b

    return _maker


@pytest.fixture
def zobov_paths_and_names():
    def _zobov_paths_and_names():
        z = ZobovElements()
        return z

    return _zobov_paths_and_names


@pytest.fixture
def svf_popcorn_paths_and_names():
    def _svf_popcorn_paths_and_names():
        s = SVFPopcornElements()
        return s

    return _svf_popcorn_paths_and_names


@pytest.fixture
def get_first_neighbor_min_distance():
    def _get_first_neighbor_min_distance(*, cloud_with_voids):
        """Calculates the minimun distance between first neighbors in a cloud.

        Parameters
        ----------
            cloud_with_voids : array
                Array of array of coordinates (x,y,z) of tracers in the cloud.
        Returns
        -------
            d_min, dmax : tuple of floats
                Minimum, maximun distances from all the distances calculated
                beween first neighbors.
        Notes
        -----
            To achieve this task, tracers (x,y,z) should be arranged based on
            closeness to (0,0,0).

            First (x,y,z) distances to (0,0,0) are calculated. Then those
            distances are sorted, ascending and then using the indexes to sort
            the tracers in the cloud, having (x0,y0,z0),...,(xN,yN,zN) sorted
            ascending.

            This ensures a value thath could be used as threshold for the
            effective_radius method where rad value is calculated
        """
        # Given a cloud calculate the distances to (0,0,0)
        dd = cdist(cloud_with_voids, np.array([[0.0, 0.0, 0.0]]))
        # Get 1D array of distances.
        dd = np.ravel(dd)
        # Get indexes that will sort the distances ascending.
        d_idx = np.argsort(dd)
        # Sort the cloud based on distance of tracers to (0,0,0)
        cwv = cloud_with_voids[d_idx]
        # Calculate the distances (now sorted) between (0,0,0) and tracers.
        nn_dist = np.linalg.norm(cwv[1:] - cwv[:-1], axis=1)
        # Get the min value from all those distances.
        d_min = np.min(nn_dist)
        # Get the max value from all those distances.
        return d_min

    return _get_first_neighbor_min_distance


@pytest.fixture
def find_bubble_neighbors():
    def _find_bubble_neighbors(*, box, cloud_with_voids, centers, rad):
        """
        Perform search of trancers inside a spherical void.

        Parameters
        ----------
            box : object
                Properties of tracers: x,y,z (coordinates) v[x,y,z] velocities
                m: total solar mass.
            cloud_with_void : nd.array
                Array of x,y,z positions of tracers.
            centers : np.array
                Array of x,y,z positions of void centers.
            rad : nd.array
                Array of void radius.
        Returns
        -------
            dist : nd.array
                Array of distances between centers and trancers.
            ind : nd.array
                Array of integer values that are the indexes of the tracers
                inside a void (regarding to the tracers index in box.)
        """
        grid = gsp.GriSPy(cloud_with_voids)
        periodic = {
            0: (box.min_, box.max_),
            1: (box.min_, box.max_),
            2: (box.min_, box.max_),
        }
        grid.set_periodicity(periodic, inplace=True)
        dist, ind = grid.bubble_neighbors(centers, distance_upper_bound=rad)
        return dist, ind

    return _find_bubble_neighbors


@pytest.fixture
def zobov_model_builder():
    def _zobov_model_builder(workdir_path=".", **kwargs):
        kwargs.setdefault("n_points", 100000)
        kwargs.setdefault("delta", -0.8)
        kwargs.setdefault(
            "centers",
            np.array([[50, 50, 50], [700, 700, 700], [500, 500, 500]]),
        )
        kwargs.setdefault("radii", 70)
        kwargs.setdefault("density_threshold", 0.1)
        kwargs_cloud = {
            key: value
            for key, value in kwargs.items()
            if key in {"seed", "lmin", "lmax", "n_points"}
        }

        cloud_with_voids_kwargs = {
            key: value
            for key, value in kwargs.items()
            if key in {"delta", "centers", "radii", "cloud"}
        }

        cloud = spherical_cloud.build_cloud(**kwargs_cloud)
        cloud_with_voids = spherical_cloud.build_spherical_void(
            cloud=cloud, **cloud_with_voids_kwargs
        )
        x, y, z = np.hsplit(cloud_with_voids, 3)
        box = Box(x=np.ravel(x), y=np.ravel(y), z=np.ravel(z))
        workdir_path = pathlib.Path(workdir_path)
        model = ZobovVF(box_size=box.max_, workdir=workdir_path)

        parameters = model.model_find(box)
        tinv, props, extra = model.build_voids(parameters)
        return extra

    return _zobov_model_builder


# =============================================================================
# POPCORN
# =============================================================================


@pytest.fixture
def find_tracers():
    def _find_tracers(properties_spheres, box):
        centers = np.array(properties_spheres[["x", "y", "z"]])
        rad = np.array(properties_spheres["radius"])

        grid = box_to_grid.get_grispy_grid_from_box(box=box)

        # For each sphere that is part of a single void, search tracers by
        # looking in sphere (bubble search)
        dist, ind = grid.bubble_neighbors(centers, distance_upper_bound=rad)

        # Map each array of tracers indexes to a set, so we don't mind
        # repeating indexes
        indx_set = np.array(list(map(set, ind)))

        # Look for tracers index belongig to spheres that are in the same void
        tracers_in_voids = [
            list(
                # 4) Perform the union of sets with tracers that belong to
                # spheres that compose a single void
                reduce(
                    set.union,
                    list(
                        # 3) Get the set indexes of tracers belongig to each
                        # sphere that compose a void
                        indx_set[
                            np.where(
                                # 2)For each unique void index filter the
                                # spheres with that index
                                properties_spheres.id
                                == idx
                            )[0]
                            # 1)Consider the unique indexes that link spheres
                            # with voids
                        ]
                    ),
                )
            )
            for idx in properties_spheres.id.unique()
        ]
        # Return list of indexes of tracers belonging to each void
        return tracers_in_voids

    return _find_tracers


@pytest.fixture
def generic_void_builder(build_box_with_eq_voids, find_bubble_neighbors):
    def _generic_void_builder(rad):
        # All voids have same radius: rad
        cloud = spherical_cloud.build_cloud()
        box, threshold, centers, cloud_with_voids = build_box_with_eq_voids(
            cloud=cloud, rad=rad
        )
        dist, idx = find_bubble_neighbors(
            box=box,
            cloud_with_voids=cloud_with_voids,
            centers=centers,
            rad=rad * np.ones(len(centers)),
        )
        vds = Voids(
            method="VoidFinderMethod",
            box=box,
            tracers_in_voids_=tuple(idx),
            centers_=centers,
            extra_={
                "radius": rad * np.ones(len(centers)),
                "properties": {
                    "VoidVol": (4 / 3)
                    * np.pi
                    * (rad * np.ones(len(centers))) ** 3
                },
            },
        )
        return vds

    return _generic_void_builder
