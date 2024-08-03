import grispy as gsp

import numpy as np

import pandas as pd

import uttr


@uttr.s(repr=False)
class VoidProperties:
    """
    A class to represent properties with various numerical attributes.

    Attributes
    ----------
    id : int
        Void Id.
    r : float
        Void Radius (32-bit).
    x : float
        Void x coordinate center (32-bit).
    y : float
        Void y coordinate center (32-bit).
    z : float
        Void z coordinate center (32-bit).
    delta_r : float
        Integrated density contrast (32-bit).

    Methods
    -------
    __repr__():
        Return a string representation of the VoidProperties object.
    """

    id = uttr.ib(converter=int)
    r = uttr.ib(converter=np.float32)
    x = uttr.ib(converter=np.float32)
    y = uttr.ib(converter=np.float32)
    z = uttr.ib(converter=np.float32)
    delta_r = uttr.ib(converter=np.float32)

    def __repr__(self):
        """
        Return a string representation of the VoidProperties object.

        Returns
        -------
        str
            String representation of the VoidProperties object.
        """
        cls_name = type(self).__name__
        return f"<{cls_name} void_number={self.id}>"


def get_void_properties(*, popcorn_output_file_path):
    """
    Read void properties from a file and return them as a tuple of
    `VoidProperties` objects.

    Parameters
    ----------
    popcorn_output_file_path : str
        Path to the file containing void properties. The file should
        contain lines where each line represents the properties of
        a void, with values separated by whitespace.

    Returns
    -------
    tuple of VoidProperties
        A tuple containing `VoidProperties` objects, one for each
        line in the file.

    Notes
    -----
    Each line in the file should contain values for the following
    attributes: 'id', 'r', 'x', 'y', 'z', 'delta_r'. The values
    should be separated by whitespace and appear in this order.

    """

    with open(popcorn_output_file_path, "r") as f:
        voids = f.readlines()
    parameters = ["id", "r", "x", "y", "z", "delta_r"]
    all_void_properties = []
    for void in voids:
        properties = dict(zip(parameters, void.split()))
        void_properties = VoidProperties(**properties)
        all_void_properties.append(void_properties)

    return tuple(all_void_properties)


def get_tracers_in_voids(*, box, popcorn_output_file_path):
    """
    Finds particles inside each void using Grispy.

    Parameters
    ----------
    box : object
        Box object (see box)

    popcorn_output_file_path : str
        Path to the file containing void properties. The file should
        have columns: 'id', 'rad', 'x', 'y', 'z', 'delta'. Each row
        represents a void with its properties including radius and
        center coordinates.

    Returns
    -------
    np.ndarray
        An array of indices of particles (refered to the Box object that are
        within each void. The shape of the array depends on the implementation
        of `grid.bubble_neighbors`.

    Notes
    -----
    The function reads the void properties from the specified file,
    including void centers and radii. It then creates a Grispy grid
    using the provided particle coordinates and finds particles that
    lie within the voids based on the radius of each void.

    The `grid.bubble_neighbors` method is used to determine which
    particles fall within the specified distance (radius) from each
    void center.
    """
    # First get x,y,z array elements and combine them into a np.array(x,y,z):
    x = box.x.value
    y = box.y.value
    z = box.z.value
    xyz = np.column_stack((x, y, z))

    # Get radius and centers
    df = pd.read_csv(
        popcorn_output_file_path,
        delim_whitespace=True,
        names=["id", "rad", "x", "y", "z", "delta"]
    )
    void_xyz = df[["x", "y", "z"]].to_numpy()
    void_rad = df["rad"].to_numpy()

    # Build grispy grid
    grid = gsp.GriSPy(xyz)

    # Get tracers
    dist, ind = grid.bubble_neighbors(void_xyz, distance_upper_bound=void_rad)
    return ind
