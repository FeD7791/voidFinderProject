import struct

import numpy as np


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
        *,
        void_properties,
        tracer_volumes,
        particle_by_voids,
        box
        ):
    """
    Calculates the center and the effective radii of each void.

    Parameters
    ----------
    void_properties : list
        List of objects VoidProperties.
    tracer_volumes : list
        List of the voronoi cell volumes holding each particle.
    particle_by_voids : list
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
    for tracers_in_void in particle_by_voids:
        center = _calculate_barycentre(
            tracers_xyz=tracers_xyz,
            tracers=tracers_in_void,
            tracer_volumes=tracer_volumes,
        )
        centers.append(center)
    return void_r_eff, centers
