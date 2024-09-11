import numpy as np

from ..utils import Bunch
from . import vsf


class Voids:
    """
    Class used to represent voids as the particles inside of them.

    Parameters
    ----------


    Attributes
    ----------
    method : str
        Name of the method used to find the voids.
    box : Box object
        Box object that holds the properties of the box.
    tracers_in_voids : tuple
        Collection of arrays that contains the IDs of the particles
        that are inside a void.
    extra : dict
        Holds extra results and information of the run.

    Methods
    -------
    method : str
        Returns the name of the method used to find the voids.
    box : Box object
        Returns a copy of the Box object containing the tracer properties.
    tracers_in_voids : tuple
        Returns the collection of arrays containing particle IDs inside voids.
    numbers_of_voids : int
        Returns the number of voids.
    extra : dict
        Returns a dictionary with extra results and information of the run.
    void_of(tracer)
        Returns indices of voids containing a specific tracer particle.

    """

    def __init__(self, *, method, box, tracers_in_voids, centers, extra):
        if len(box) <= len(tracers_in_voids):
            raise ValueError(
                "Number of box must be lesser than the numbers of voids"
            )

        self._method = str(method)
        self._tracers = box.copy()  # the box
        self._tracers_in_voids = tuple(tracers_in_voids)  # tuple of arrays
        self._centers = np.array(centers)
        self._extra = Bunch("extra", extra)  # dict with zaraza

    @property
    def method(self):
        """str: Name of the method used to find the voids."""
        return self._method

    @property
    def box(self):
        """Box object: Box object that holds the properties of the box."""
        return self._tracers

    @property
    def centers(self):
        return self._centers

    @property
    def tracers_in_voids_(self):
        """tuple: Collection of arrays that contains the IDs of particles
        inside voids."""
        return self._tracers_in_voids

    @property
    def numbers_of_voids_(self):
        """int: Number of voids."""
        return len(self._tracers_in_voids)

    @property
    def extra_(self):
        """dict: Holds extra results and information of the run."""
        return dict(self._extra)

    e_ = extra_

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
        for idx, void in enumerate(self._tracers_in_voids):
            if tracer in void:
                voids_w_tracer.append(idx)
        return np.array(voids_w_tracer)

    def effective_radius(self, *, delta=-0.9, n_neighbors=100, n_cells=64):
        return vsf.effective_radius(
            self.centers,
            self.box,
            delta=delta,
            n_neighbors=n_neighbors,
            n_cells=n_cells,
        )
