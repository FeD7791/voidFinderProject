from collections.abc import Sequence
import dataclasses as dtclss
import warnings

import grispy as gsp

import numpy as np


class EffectiveRadiusErrors:
    """Enum of effective radius calculation errors."""

    NO_ERROR = 0
    MAYBE_NEAR_ANOTHER_VOID = 1
    EXEED_CRITICAL = 2
    UNDER_CRITICAL = 3


@dtclss.dataclass(frozen=True, slots=True, repr=False)
class _EffectiveRadius(Sequence):
    """
    A dataclass representing the effective radius calculation results.

    Parameters
    ----------
    delta : float
        The density contrast parameter.
    n_neighbors : int
        The number of neighbors considered in the calculation.
    n_cells : int
        The number of cells used in the spatial grid.
    errors : np.ndarray
        Array of error codes for each calculation.
    radius : np.ndarray
        Array of calculated effective radii.
    tracers : np.ndarray
        Array of tracer particles for each void.
    densities : np.ndarray
        Array of density values for each void.

    """

    delta: float
    n_neighbors: int
    n_cells: int
    errors: np.ndarray
    radius: np.ndarray
    tracers: np.ndarray
    densities: np.ndarray

    @property
    def argerrors(self):
        """Return a boolean array indicating which calculations resulted in \
        errors.

        Returns
        -------
        np.ndarray
            Boolean array where True indicates an error occurred.
        """
        return self.errors != EffectiveRadiusErrors.NO_ERROR

    def __repr__(self):
        """Returns a string representation of the dataclass.

        Returns
        -------
        str
            The string representation of the dataclass.
        """
        delta = self.delta
        n_neighbors = self.n_neighbors
        n_cells = self.n_cells
        good = np.sum(~self.argerrors)
        total = len(self)
        return (
            "<effective_radius "
            f"{delta=} {n_neighbors=} {n_cells=} | {good}/{total}>"
        )

    def __len__(self):
        """Returns the length of the dataclass.

        Returns
        -------
        int
            The length of the dataclass.
        """
        return len(self.errors)

    def __getitem__(self, slicer):
        """
        Returns a subset of the dataclass.

        Parameters
        ----------
        slicer : slice
            The slice to apply to the dataclass.

        Returns
        -------
        tuple
            The subset of the dataclass.

        """
        return (
            self.errors.__getitem__(slicer),
            self.radius.__getitem__(slicer),
            self.tracers.__getitem__(slicer),
            self.densities.__getitem__(slicer),
        )


def _sigle_void_eradius(idx, n_neighbors, crit_density, distance, nn):
    """Calculate the effective radius for a single void.

    Parameters
    ----------
    idx : int
        Index of the void center.
    n_neighbors : int
        Number of neighbors to consider.
    crit_density : float
        Critical density threshold.
    distance : np.ndarray
        Array of distances to neighboring particles.
    nn : np.ndarray
        Array of nearest neighbor indices.

    Returns
    -------
    tuple
        A tuple containing:
        (error_code, effective_radius, void_tracers, void_density).

    """
    # Find density values for n_nat particles at radius d
    n_nat = np.arange(1, n_neighbors + 1)
    density_n_nat_d = (3 * n_nat) / (4 * np.pi * distance**3)

    # Find all density values that are less than crit_density
    dens_values = np.where(density_n_nat_d < crit_density)[0]

    # This means that all calculated densities are above crit density,
    # probably not a void
    if len(dens_values) == 0:
        # void_error, void_radius, void_tracers, void_density
        return (
            EffectiveRadiusErrors.EXEED_CRITICAL,
            np.nan,
            [],
            density_n_nat_d,
        )

    elif len(dens_values) == len(density_n_nat_d):
        # warning
        warnings.warn(
            f"All values under critical Density for center {idx}",
            RuntimeWarning,
        )
        # void_error, void_radius, void_tracers, void_density
        return (
            EffectiveRadiusErrors.UNDER_CRITICAL,
            np.nan,
            [],
            density_n_nat_d,
        )

    else:
        # From the values that fulfill the latter condition find the index
        # of the value with max radii
        dist_max_index = np.where(distance == np.max(distance[dens_values]))[
            0
        ][0]

        # This would mean that the biggest
        # radius is asociated to a density
        # that has yet not crosed the crit density You either have to
        # increase n_neighbors or assume this is a 'noisy' void
        if distance[dist_max_index] == distance[-1]:
            # void_error, void_radius, void_tracers, void_density
            return (
                EffectiveRadiusErrors.MAYBE_NEAR_ANOTHER_VOID,
                np.nan,
                [],
                density_n_nat_d,
            )

        # Final radii is half distance between distk_max_index and
        # dist_max_index +1
        else:
            # void_error, void_radius, void_tracers, void_density
            radius = (
                distance[dist_max_index + 1] + distance[dist_max_index]
            ) / 2
            tracers = nn[idx][:dist_max_index]

            # void_error, void_radius, void_tracers, void_density
            return (
                EffectiveRadiusErrors.NO_ERROR,
                radius,
                tracers,
                density_n_nat_d,
            )


def effective_radius(centers, box, *, delta, n_neighbors, n_cells):
    """Calculate the effective radius for multiple void centers.

    Parameters
    ----------
    centers : array-like
        Array of void center coordinates.
    box : object
        Object representing the simulation box, with x, y, z attributes and a
        size() method.
    delta : float, optional
        Density contrast parameter.
    n_neighbors : int, optional
        Number of neighbors to consider.
    n_cells : int, optional
        Number of cells for spatial gridding.

    Returns
    -------
    effectiveradius
        An object containing the results of the effective radius calculations.

    """
    # Create spatial gridding
    xyz = np.column_stack((box.x.value, box.y.value, box.z.value))
    grid = gsp.GriSPy(xyz, copy_data=False, N_cells=n_cells)

    # For each center, get the distance for the n nearest tracers and their
    # index
    distances, nn = grid.nearest_neighbors(centres=centers, n=n_neighbors)

    tracers = np.zeros(len(distances), dtype=object)
    radius = np.zeros(len(distances), dtype=float)
    errors = np.zeros(len(distances), dtype=int)
    densities = np.zeros(len(distances), dtype=object)

    # This is the density below which all voids should be
    # to be considered an underdensity
    crit_density = (1 + delta) * (len(box) / (box.size() ** 3))

    # Find the effective radius for each center
    for idx, distance in enumerate(distances):
        verror, vradius, vtracers, vdensity = _sigle_void_eradius(
            idx=idx,
            n_neighbors=n_neighbors,
            crit_density=crit_density,
            distance=distance,
            nn=nn,
        )

        errors[idx] = verror
        radius[idx] = vradius
        tracers[idx] = vtracers
        densities[idx] = vdensity

    # create effective radius object
    eradius = _EffectiveRadius(
        delta=delta,
        n_neighbors=n_neighbors,
        n_cells=n_cells,
        errors=errors,
        radius=radius,
        tracers=tracers,
        densities=densities,
    )

    return eradius


# =============================================================================
# VSF
# =============================================================================
