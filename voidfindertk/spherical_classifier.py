import scipy
import numpy as np
from attrs import define
from .box import Box
from .spherical_voids import SphericalVoids


@define(repr=False)
class Box_SphericalVoids_Classifier(Box, SphericalVoids):
    """Inherits from both Box and ShpericalVoids classes and provides a method
    to create a sparse matrix representing the relationship between spherical
    voids and tracer particles within a simulation box.

    Methods
    -------
    _sparse_matrix(self, tol=0.0)
        Generates a sparse CSR matrix indicating which voids are potentially
        influencing each tracer particle based on their relative positions.

    Parameters
    ----------
    self : Box_SphericalVoids_Classifier object
        The instance of the class.
    tol : float, optional
        Tolerance level for considering a void to be influencing a tracer
        particle. Defaults to 0.0.

    Returns
    -------
    scipy.sparse.csr_matrix
        A sparse CSR matrix where each row represents a tracer particle and
        each column represents a void. The value at a specific row-column
        intersection is 1.0 if the corresponding void is potentially
        influencing the tracer particle (Tracer inside the void) based on
        the distance between their centers and the void's radius,
        otherwise 0.0.
    """

    def _sparse_matrix(self, tol=0.0):
        tolerance = tol * np.ones(self.rad.shape)
        rad = np.array(self.rad) + tolerance
        pos_void = np.array([self.x_void, self.y_void, self.z_void])
        pos_void = np.transpose(pos_void)
        arr = []
        for i in range(self._len):
            pos_box = np.array(
                [[self.x[i].value, self.y[i].value, self.z[i].value]]
            )
            d = scipy.spatial.distance.cdist(
                pos_box, pos_void, metric="euclidean"
            )  # Calculates the distance from a tracer to each void center
            rad_center_distance_comparing = np.greater(rad, d[0]).astype(
                float
            )  # Compares radius of a void and distance from the void center to the particle 1 if rad>d # noqa: E501
            arr.append(rad_center_distance_comparing)
        # Transform to sparse matrix
        output_sparse_matrix = scipy.sparse.csr_matrix(arr)
        return output_sparse_matrix

    def __repr__(self):
        """Representation method.

        Returns
        -------
            str
                Size of Box and calculated Voids
        """
        cls_name = type(self).__name__
        # length_box = self._len
        # length_voids = self._void_len
        return f"<{cls_name} >"
