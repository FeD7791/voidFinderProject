import dataclasses as dtclss
import warnings
from collections.abc import Sequence

import grispy as gsp

import numpy as np


class EffectiveRadiusErrors:
    """
    Enumeration for error codes related to the computation of effective radii
    of voids.

    This class defines a set of error codes used to indicate the status or
    issues encountered during the calculation of effective radii for voids in
    a system. These codes help in diagnosing and understanding the results of
    the radius computation process.

    Attributes
    ----------
    NO_ERROR : int
        Error code indicating that the effective radius computation completed
        successfully without issues.
    MAYBE_NEAR_ANOTHER_VOID : int
        Error code indicating that the void might be near another void,
        leading to potential inaccuracies.
    EXEED_CRITICAL : int
        Error code indicating that the computed densities exceeded the
        critical density, suggesting the region is not considered a void.
    UNDER_CRITICAL : int
        Error code indicating that all computed densities are below the
        critical density, which imply that the number of nearest neighbors
        parameter should be increased to find the right radius of the void.

    """

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
    """
    Compute the effective radius of voids based on nearest neighbor distances
    using grispy.

    This function calculates the effective radius of voids by considering
    the distances to the nearest neighbors for each center and comparing
    against a critical density. It returns the effective radius, error,
    tracer particles within the void, and density map for each void.

    Parameters
    ----------
    centers : array-like, shape (n_centers, 3)
        Coordinates of the centers for which the effective radius is to be
        computed.
    box : Box
        Box object containing the properties of the spatial domain.
    delta : float, optional
        Parameter to adjust the critical density. The critical density is
        calculated as (1 + delta) * (number of tracers / (box volume)^3).
        Default is -0.9.
    n_neighbors : int, optional
        Number of nearest neighbors to consider for each center.
        Default is 100.
    n_cells : int, optional
        Number of cells used for the spatial grid. Default is 64.

    Returns
    -------
    errors : list of float
        List of errors associated with each void, posble values are:
        0 : No error
        1 : Local overdensity, maybe related to two near underdensities.
        2 : Densitie map over the crital density minima. Probably not a void
        3 : Densitie map under the critical density minima. Increase the
        number of nearest neighbors used to perform the search.
    radius : list of float
        List of effective radii for each void.
    tracers : tuple of numpy.ndarray
        Tuple where each element is an array containing the IDs of the tracers
        within the corresponding void.
    densities : list of float
        List of densities for each void.

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
