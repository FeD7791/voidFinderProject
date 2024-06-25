import numpy as np


class Voids:
    """
    A class representing voids and their associated tracers.

    Parameters
    ----------
    tracers : iterable
        List or array-like object containing tracers associated with voids.
    voids : iterable
        List or array-like object containing voids.

    Attributes
    ----------
    tracers : Object
        Box object that contains tracers and their properties.
    voids : list
        List of ZobovVoids objects.

    Methods
    -------
    voids_numbers()
        Returns the number of voids.
    void_of(tracer)
        Finds and returns indices of voids containing a specific tracer.

    Notes
    -----
    This class inherits from ModelVoid.

    """

    def __init__(self, *, method, tracers, voids, extra):
        if len(tracers) <= len(voids):
            raise ValueError()

        self._method = str(method)
        self._tracers = tracers  # the tracers
        self._voids = tuple(voids)  # tuple of arrays
        self._extra = dict(extra)  # dict with zaraza

    @property
    def method(self):
        return self._method

    @property
    def tracers(self):
        return self._tracers

    @property
    def voids(self):
        return self._voids

    @property
    def numbers_of_voids(self):
        return len(self.voids)

    @property
    def extra(self):
        return dict(self._extra)

    # REPR ====================================================================
    def __repr__(self):
        return ""

    # utilities ===============================================================

    def void_of(self, tracer):
        voids_w_tracer = []
        for idx, void in enumerate(self.voids):
            if tracer in void.tracers_in_void:
                voids_w_tracer.append(idx)
        return np.array(voids_w_tracer)
