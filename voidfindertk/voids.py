import numpy as np

from .utils import Bunch


class Voids:
    """
    Class used to represent voids as the particles inside of them.

    Parameters
    ----------
        method : str
            Name of the Method used to find the voids.
        tracers : Box object
            Box object that holds the properties of the tracers.
        voids : tuple
            Collection of arrays that contains the IDs of the particles
            that are inside a void. Each element holds the IDs of particles in
            a void in ascending void number (see VoidProperties
            class in _postprocessing module).
        extra :
            Holds extra results and information of the run.

    Attributes
    ----------
    method : str
        Name of the method used to find the voids.
    tracers : Box object
        Box object that holds the properties of the tracers.
    voids : tuple
        Collection of arrays that contains the IDs of the particles
        that are inside a void.
    extra : dict
        Holds extra results and information of the run.

    Methods
    -------
    method : str
        Returns the name of the method used to find the voids.
    tracers : Box object
        Returns a copy of the Box object containing the tracer properties.
    voids : tuple
        Returns the collection of arrays containing particle IDs inside voids.
    numbers_of_voids : int
        Returns the number of voids.
    extra : dict
        Returns a dictionary with extra results and information of the run.
    void_of(tracer)
        Returns indices of voids containing a specific tracer particle.

    """

    def __init__(self, *, method, tracers, voids, extra):
        if len(tracers) <= len(voids):
            raise ValueError(
                "Number of tracers must be lesser than the numbers of voids"
            )

        self._method = str(method)
        self._tracers = tracers.copy()  # the tracers
        self._voids = tuple(voids)  # tuple of arrays
        self._extra = Bunch("extra", extra)  # dict with zaraza

    @property
    def method(self):
        """str: Name of the method used to find the voids."""
        return self._method

    @property
    def tracers(self):
        """Box object: Box object that holds the properties of the tracers."""
        return self._tracers

    @property
    def voids_(self):
        """tuple: Collection of arrays that contains the IDs of particles
        inside voids."""
        return self._voids

    @property
    def numbers_of_voids_(self):
        """int: Number of voids."""
        return len(self._voids)

    @property
    def extra_(self):
        """dict: Holds extra results and information of the run."""
        return dict(self._extra)

    e_ = extra_

    # REPR ====================================================================
    def __repr__(self):
        return (
            f"<Voids '{self.method}' "
            f"{self.numbers_of_voids_}V, {len(self.tracers)}T>"
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
        for idx, void in enumerate(self._voids):
            if tracer in void:
                voids_w_tracer.append(idx)
        return np.array(voids_w_tracer)
