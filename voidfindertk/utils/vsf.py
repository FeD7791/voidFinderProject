import warnings

import grispy as gsp

import numpy as np


class EffectivRadiusErrors:
    NO_ERROR = 0
    MAYBE_NEAR_ANOTHER_VOID = 1
    EXEED_CRITICAL = 2
    UNDER_CRITICAL = 3


def _void_effr(idx, n_neighbors, crit_density, distance, nn):

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
            EffectivRadiusErrors.EXEED_CRITICAL,
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
            EffectivRadiusErrors.UNDER_CRITICAL,
            np.nan,
            [],
            density_n_nat_d,
        )

    else:
        # From the values that fulfill the latter condition find the index
        # of the value with max radii
        dist_max_index = np.where(distance == np.max(distance[dens_values]))[0][0]

        # This would mean that the biggest
        # radius is asociated to a density
        # that has yet not crosed the crit density You either have to
        # increase n_neighbors or assume this is a 'noisy' void
        if distance[dist_max_index] == distance[-1]:
            # void_error, void_radius, void_tracers, void_density
            return (
                EffectivRadiusErrors.MAYBE_NEAR_ANOTHER_VOID,
                np.nan,
                [],
                density_n_nat_d,
            )

        # Final radii is half distance between distk_max_index and
        # dist_max_index +1
        else:
            # void_error, void_radius, void_tracers, void_density
            radius = (distance[dist_max_index + 1] + distance[dist_max_index]) / 2
            tracers = nn[idx][:dist_max_index]

            # void_error, void_radius, void_tracers, void_density
            return (
                EffectivRadiusErrors.NO_ERROR,
                radius,
                tracers,
                density_n_nat_d,
            )


def effective_radius(centers, box, *, delta=-0.9, n_neighbors=100, n_cells=32):
    """
    Calculates the radius for a void when the center and box of tracers are
    provided.

    The final radius is the arithmetic mean between density(n,r_n) and
    density(n+1,r_n+1), where density(n,r_n) < crit_density < density(n+1,
    r_n+1)
    where density(n,r_n) = n/(4*pi*r_n**3)/3 (Spherical Volume)

    Parameters
    ----------
        centers : array of (x,y,z)
        Array of tracers coordinates.

        n_neighbors : int
        Maximun number of tracers in each void used to perform the search.

        box : Box Object
        Object with the tracer properties.

        delta : float -1< delta < 0
        Integrated density contrast.

    Returns
    -------
        rad_new : numpy array
        Radius calculated for each void.The radius could have 4 values:
            0 : Radius so thath density < crit_density, but had some radius
                that had density > crit_density.
            -2 : density < crit_density it is necesary to increment number
                of neighbors
            -1 : All radii have that crit_density < density. Probably not
                a void.

        tracers_new : tuple of lists
        List of indexes listing the tracers inside each void

        density_values : list of arrays
        For each void gives 1 < n <= n_neighbors density(n,r_n) calculations.
    """

    xyz = np.column_stack((box.x.value, box.y.value, box.z.value))
    grid = gsp.GriSPy(xyz, copy_data=False, n_cells=n_cells)

    # For each center, get the distance for the n nearest tracers and their
    # index
    distances, nn = grid.nearest_neighbors(centres=centers, n=n_neighbors)

    tracers = []
    radius = []
    errors = []
    densities = []

    # es la densidad por la cual todos los voids deberian estar por debajo
    # para considerarse una subdensidad
    crit_density = (1 + delta) * (len(box) / (box.size() ** 3))

    for idx, distance in enumerate(distances):
        void_error, void_radius, void_tracers, void_density = _void_effr(
            idx=idx,
            n_neighbors=n_neighbors,
            crit_density=crit_density,
            distance=distance,
            nn=nn,
        )

        errors.append(void_error)
        radius.append(void_radius)
        tracers.append(void_tracers)
        densities.append(void_density)

    return errors, radius, tuple(tracers), densities


# =============================================================================
# VSF
# =============================================================================

