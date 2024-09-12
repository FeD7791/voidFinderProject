import numpy as np

import attrs

from ..utils import Bunch
from . import vsf, plot_acc
from .box import Box


@attrs.define(frozen=True, repr=False)
class Voids:
    """
    A class to represent and manage voids in a system of particles.

    This class provides methods to analyze voids and the particles within them,
    including retrieving the method used to find the voids, accessing the box
    properties, and finding specific voids associated with tracer particles.

    Parameters
    ----------
    method : str
        The name of the method used to find the voids.
    box : Box
        The Box object containing properties of the box.
    tracers_in_voids : tuple of numpy.ndarray
        A tuple of arrays, where each array contains the IDs of particles
        inside a void.
    centers : array-like
        Coordinates of the centers of the voids.
    extra : dict
        Additional results and information of the run.

    Attributes
    ----------
    method : str
        The name of the method used to find the voids.
    box : Box
        A copy of the Box object containing tracer properties.
    centers : numpy.ndarray
        Array of coordinates of the centers of the voids.
    tracers_in_voids : tuple of numpy.ndarray
        Collection of arrays that contains the IDs of particles inside voids.
    numbers_of_voids : int
        The number of voids.
    extra : dict
        Dictionary with additional results and information of the run.

    Methods
    -------
    void_of(tracer) :
        Returns indices of voids containing a specific tracer particle.
    all_effective_radius(delta=-0.9, n_neighbors=100, n_cells=64) :
        Computes the effective radius for all voids.
    effective_radius(void_idx, delta=-0.9, n_neighbors=100, n_cells=64) :
        Computes the effective radius of a specific void.

    Raises
    ------
    ValueError
        If the number of boxes is not lesser than the number of voids.

    Notes
    -----
    - The `extra` attribute is managed using a `Bunch` object, which is
    essentially a dictionary.
    - The `tracers_in_voids` attribute is a tuple of arrays where each array
    contains IDs of particles inside a void.

    """

    # came from finde and data
    method: str = attrs.field(converter=str)
    box: Box = attrs.field(converter=Box.copy)

    # if end with "_" is calculated
    tracers_in_voids_: tuple = attrs.field(converter=tuple)
    centers_: np.ndarray = attrs.field(converter=np.array)
    extra_: Bunch = attrs.field(converter=lambda e: Bunch("extra", e))

    # plot accessor
    plot: plot_acc.VoidPlotter = attrs.field(
        init=False,
        default=attrs.Factory(plot_acc.VoidPlotter, takes_self=True),
    )

    def __attrs_post_init__(self):
        if len(self.box) <= len(self.tracers_in_voids_):
            raise ValueError(
                "Number of box must be lesser than the numbers of voids"
            )

    @property
    def numbers_of_voids_(self):
        """int: Number of voids."""
        return len(self.tracers_in_voids_)

    @property
    def e_(self):
        """dict: Holds extra results and information of the run."""
        return self.extra_

    # REPR ====================================================================
    def __repr__(self):
        return (
            f"<Voids '{self.method}' "
            f"{self.numbers_of_voids_}V, {len(self.box)}T>"
        )

    # utilities ===============================================================

    def void_of(self, tracer):
        """
        Returns indices of voids containing a specific tracer particle.

        Parameters
        ----------
        tracer : int
            ID of the tracer particle to search for in voids.

        Returns
        -------
        numpy.ndarray
            Array of indices of voids containing the tracer particle.

        """
        voids_w_tracer = []
        for idx, void in enumerate(self.tracers_in_voids):
            if tracer in void:
                voids_w_tracer.append(idx)
        return np.array(voids_w_tracer)

    def effective_radius(self, *, delta=-0.9, n_neighbors=100, n_cells=64):
        """Calculate the effective radius of the voids.

        Parameters
        ----------
        delta : float
            The delta parameter for the effective radius calculation.
            Defaults to -0.9.
        n_neighbors : int
            The number of neighbors to consider for the effective radius
            calculation. Defaults to 100.
        n_cells : int
            The number of cells to consider for the effective radius
            calculation. Defaults to 64.

        Returns
        -------
        effective_radius :
            A dataclass containing the effective radius results.

        """

        return vsf.effective_radius(
            self.centers_,
            self.box,
            delta=delta,
            n_neighbors=n_neighbors,
            n_cells=n_cells,
        )

    def void_size_function(
        self,
        *,
        scale_1_num_samples=7,
        scale_2_num_samples=2,
        **kwargs,
    ):
        # Get the radii
        effective_radius = self.effective_radius(**kwargs)

        log_of_radius, count, delta = vsf.void_size_function(
            effective_radius=effective_radius,
            box=self.box,
            scale_1_num_samples=scale_2_num_samples,
            scale_2_num_samples=scale_2_num_samples,
        )

        return log_of_radius, count, delta

    vsf = void_size_function
