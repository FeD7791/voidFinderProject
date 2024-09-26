import struct

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
        Id of void
    x : float
        Void x coordinate center (32-bit).
    y : float
        Void y coordinate center (32-bit).
    z : float
        Void z coordinate center (32-bit).
    r : float
        Void Radius (32-bit).
    density_contrast : float
        Spherically-averaged density contrast (32-bit).

    Methods
    -------
    __repr__():
        Return a string representation of the VoidProperties object.
    """

    id = uttr.ib(converter=int)
    x = uttr.ib(converter=np.float32)
    y = uttr.ib(converter=np.float32)
    z = uttr.ib(converter=np.float32)
    r = uttr.ib(converter=np.float32)
    density_contrast = uttr.ib(converter=np.float32)

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


def read_volume_file(*, filename):
    """
    Read volume data from the binary volume file and return it as a NumPy
    array.

    Parameters
    ----------
    filename : str
        The path to the binary volume file containing the volume data. This
        file contains the volumes of the voronoi cells where each particle is
        in.

    Returns
    -------
    volumes : numpy.ndarray
        A 1-D NumPy array of type `np.float32` containing the volume data read
        from the binary file.

    Notes
    -----
    The indexes of the returned array are directly related to the tracers
    index in the box object.

    Examples
    --------
    >>> volumes = read_volume_file(filename='voloutput_vozinit.dat')
    >>> print(volumes)
    [1.234 5.678 9.101]
    """
    with open(filename, "rb") as f:
        # Read number of tracers
        number_voids = struct.unpack("i", f.read(4))[0]
        # Read volumes
        volumes = np.zeros(number_voids, dtype=np.float32)
        for i in range(number_voids):
            volume = struct.unpack("d", f.read(8))[0]
            volumes[i] = np.float32(volume)
    return volumes


def _get_tracers_xyz(*, box):
    """
    Gets the x,y,z coordinates of the tracers inside the box.

    Parameters
    ----------
    box : Object Box
        Object that holds the properties of the input tracers

    Returns
    -------
    xyz_arr : numpy.array
        Array of N rows, 3 cols where each row is the x,y,z position of a
        tracer.
    """
    tracer_x = box.x.value
    tracer_y = box.y.value
    tracer_z = box.z.value
    xyz_arr = np.stack([tracer_x, tracer_y, tracer_z], axis=1)
    return xyz_arr


def _calculate_barycentre(*, tracers_xyz, tracers, tracer_volumes):
    """
    Calucates the barycentre of a single void. A void is composed of several
    voronoi cells. Each particle is inside a unique voronoi cell.

    The barycentre of a particular void iscalculated as the weigthed sum of
    the particles position times the volume of its voronoi cell divided the
    total volume of the void.

    Parameters
    ----------
    tracers_xyz : numpy.array
        Array of [x,y,z] tracers positions
    tracers : list
        List of indexes of each tracer that belongs to a particular void.
    tracer_volumes : numpy.array
        Array [v1,v2,v3...] of voronoi cell volumes of each particle.

    Returns
    -------
    center : numpy.array
        [X,Y,Z] coordinates of the center of the void.


    """
    tracer_volumes = tracer_volumes[tracers]
    arr = tracers_xyz[tracers]
    center = np.average(arr, weights=tracer_volumes, axis=0)
    return center


def _get_volumes_from_properties(*, void_properties):
    """
    Get the total volume of all the voids in the properties input.

    Parameters
    ----------
    void_properties : list of Objects VoidProperties
        List of objects that have the properties of the ZOBOV voids
        properties. Concretely the void_vol property is extracted from all
        the list.

    Returns
    -------
    void_volumes : list
        List of void_vol property of each void.
    """
    void_volumes = []
    for void in void_properties:
        void_volumes.append(void.void_vol.value)
    void_volumes = np.array(void_volumes, dtype=np.float32)
    return void_volumes


def _calculate_r_eff(*, void_volumes):
    """
    Get the effective Radii for a given volume

    Parameters
    ----------
    void_volumes : value/array
    A single volume value or array of volume values.

    Returns
    -------
    r_eff : value/array
    Effective radius for each provided volume.
    """
    r_eff = ((3 / (4 * np.pi)) * void_volumes) ** (1 / 3)
    return r_eff


def get_center_and_radii(
    *, void_properties, tracer_volumes, tracers_in_voids, box
):
    """
    Calculates the center and the effective radii of each void.

    Parameters
    ----------
    void_properties : list
        List of objects VoidProperties.
    tracer_volumes : list
        List of the voronoi cell volumes holding each particle.
    tracers_in_voids : list
        List of particles within voids.
    box : Object
        Box Object that holds information about the tracers data.

    Returns
    -------
    void_r_eff : array
        Array of effective radius of each void.
    centers : list
        List of (x,y,z) coordinates of each void center.

    """
    # Get the void_volumes
    void_volumes = _get_volumes_from_properties(
        void_properties=void_properties
    )
    # Get r_eff
    void_r_eff = _calculate_r_eff(void_volumes=void_volumes)

    # Get tracers xyz coords
    tracers_xyz = _get_tracers_xyz(box=box)

    # Get centers
    centers = []
    for tracers_in_void in tracers_in_voids:
        center = _calculate_barycentre(
            tracers_xyz=tracers_xyz,
            tracers=tracers_in_void,
            tracer_volumes=tracer_volumes,
        )
        centers.append(center)
    centers = np.array(centers)
    return void_r_eff, centers


def save_r_eff_center(*, centers, r_eff, path):
    """
    Save centers and effective radii to a tab-separated values file.

    Parameters
    ----------
    centers : array-like, shape (n_samples, 3)
        An array or list of coordinates representing the centers. Each entry
        should be a sequence of three values, corresponding to x, y, and z
        coordinates.

    r_eff : array-like, shape (n_samples,)
        An array or list of effective radii, where each value corresponds to
        the effective radius for the respective center.

    path : str
        The file path where the tab-separated values file will be saved.

    Returns
    -------
    None
        This function does not return any value.

    Notes
    -----
    The output file will have columns for x, y, z coordinates and r_eff,
    separated by tabs. The file will not include an index or header row.
    """
    df = pd.DataFrame(centers)
    df.columns = ["x", "y", "z"]
    df["r_eff"] = r_eff
    df.to_csv(path, index=False, header=False, sep="\t")


def save_xyz_tracers(*, box, path):
    """
    Save x, y, and z coordinates from a box object to a tab-separated values
    file.

    Parameters
    ----------
    box : object
        An object with attributes `x`, `y`, and `z`, each of which should have
        a `value` attribute that is an array-like sequence of coordinates.

    path : str
        The file path where the tab-separated values file will be saved.

    Returns
    -------
    None
        This function does not return any value.

    Notes
    -----
    The output file will have columns for x, y, and z coordinates, separated
    by tabs. The file will not include an index or header row.
    """
    x = box.x.value
    y = box.y.value
    z = box.z.value
    xyz = np.column_stack((x, y, z))
    df = pd.DataFrame(xyz)
    df.to_csv(path, index=False, header=False, sep="\t")


def get_tracers_and_center(*, box, cbl_cleaned_path):
    """
    Finds particles inside each void using Grispy.

    Parameters
    ----------
    box : object
        Box object (see box)

    cbl_cleaned_path : str
        Path to the file containing void properties. The file should
        have columns: 'x', 'y', 'z', 'rad','delta'. Each row
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
        cbl_cleaned_path,
        delim_whitespace=True,
        names=["x", "y", "z", "rad", "delta"],
    )
    void_xyz = df[["x", "y", "z"]].to_numpy()
    void_rad = df["rad"].to_numpy()

    # Build grispy grid
    grid = gsp.GriSPy(xyz)

    # Get tracers
    dist, ind = grid.bubble_neighbors(void_xyz, distance_upper_bound=void_rad)
    return ind, void_xyz


def get_dive_void_properties(*, cleaned_catalogue_path):
    with open(cleaned_catalogue_path, "r") as f:
        voids = f.readlines()
    parameters = ["id", "x", "y", "z", "r", "density_contrast"]
    all_void_properties = []
    index = 0
    for void in voids:
        void = void.split()
        void.insert(0, index)  # Add index
        properties = dict(zip(parameters, np.array(void, dtype=np.float32)))
        all_void_properties.append(VoidProperties(**properties))
        index = index + 1
    return tuple(all_void_properties)
