#!/usr/bin/env python3
# =============================================================================
# -*- coding: utf-8 -*-
# Copyright (c) 2023, Bustillos Federico, Gualpa Sebastian, Cabral Juan,
# Paz Dante, Ruiz Andres, Correa Carlos
# License: MIT
# Full Text: https://github.com/FeD7791/voidFinderProject/blob/dev/LICENSE.txt
# All rights reserved.
# =============================================================================
"""Contains functions to parse output files from the ZOBOV void finder."""

import ctypes

from astropy import units as u

import numpy as np

import uttr


@uttr.s(repr=False)
class VoidProperties:
    """
    Properties of voids found with ZOBOV.

    This class represents void properties from ZOBOV's ASCII output file.

    Attributes
    ----------
    void_number : int
        Void rank in decreasing order of VoidDensContrast.
    file_void_number : int
        Void number (starting from 0) in the first two files.
    core_particle : int
        Particle number (from 0) of the void's (and zone's) core particle.
    core_dens : float
        Core particle density, in units of mean density.
    zone_vol : float
        Central zone volume, in units of mean-density particle volume.
    zone_number_part : int
        Number of particles in the void's central zone.
    void_number_zones : int
        Number of zones in the void.
    void_vol : float
        Void volume, in units of mean-density particle volume.
    void_number_part : int
        Number of particles in the void.
    void_dens_contrast : float
        Void density contrast (ratio of critical flow density to min density).
    void_prob : float
        Probability of DensContrast arising from Poisson noise (ZOBOV eq. 1).

    Notes
    -----
    void_prob is based on DensContrast distribution from Poisson particles.

    """

    void_number = uttr.ib(converter=int)
    file_void_number = uttr.ib(converter=int)
    core_particle = uttr.ib(converter=int)
    core_dens = uttr.ib(converter=np.float32)
    zone_vol = uttr.ib(converter=np.float32, unit=u.Mpc**3)
    zone_number_part = uttr.ib(converter=int)
    void_number_zones = uttr.ib(converter=int)
    void_vol = uttr.ib(converter=np.float32, unit=u.Mpc**3)
    void_number_part = uttr.ib(converter=int)
    void_dens_contrast = uttr.ib(converter=np.float32)
    void_prob = uttr.ib(converter=np.float32)

    def __repr__(self):
        """
        Return a string representation of the VoidProperties object.

        Returns
        -------
        str
            String representation of the VoidProperties object.
        """
        cls_name = type(self).__name__
        return f"<{cls_name} void_number={self.void_number}>"


# =============================================================================
# FUNCTIONS
# =============================================================================
def parse_zones_in_void_output(
    *, executable_path, input_file_path, output_file_path
):
    """
    Parse tracers in zones output using a C library.

    Parameters
    ----------
    executable_path : str
        Path to the C library executable.
    input_file_path : str
        Path to the input file with tracers in zones data.
    output_file_path : str
        Path where the parsed output will be saved.

    Notes
    -----
    Uses ctypes to call a C function for parsing tracers in zones output.
    """
    # Get library
    clibrary = ctypes.CDLL(str(executable_path), mode=ctypes.RTLD_GLOBAL)

    # Input argtypes
    clibrary.process_files.argtypes = [ctypes.c_char_p, ctypes.c_char_p]

    # Call function
    clibrary.process_files(
        str(input_file_path).encode(), str(output_file_path).encode()
    )


def parse_tracers_in_zones_output(
    *, executable_path, input_file_path, output_file_path
):
    """
    Parse tracers in zones output using a C library.

    Parameters
    ----------
    executable_path : str
        Path to the C library executable.
    input_file_path : str
        Path to the input file with tracers in zones data.
    output_file_path : str
        Path where the parsed output will be saved.

    Notes
    -----
    Uses ctypes to call a C function for parsing tracers in zones output.
    """
    # Get library
    clibrary = ctypes.CDLL(str(executable_path), mode=ctypes.RTLD_GLOBAL)

    # Input argtypes
    clibrary.get_tracers_in_zones.argtypes = [ctypes.c_char_p, ctypes.c_char_p]

    # Call function
    clibrary.get_tracers_in_zones(
        str(input_file_path).encode(), str(output_file_path).encode()
    )


def _get_particles_in_zones(*, particles_in_zones_path):
    """
    Extract particles in voids from a parsed file.

    Parameters
    ----------
    particles_in_zones_path : str
        Path to the file with parsed particles in zones data.

    Returns
    -------
    dict
        Dictionary with zone IDs as keys and numpy arrays of particle IDs
        as values.
    """
    with open(particles_in_zones_path, "r") as f:  # Read Parsed file
        zones_particles = f.readlines()

    particles_in_zones = {}  # Dict that will contain the list of particles
    for i in np.arange(len(zones_particles)):
        # Deal with the format of the particles in zones file
        if zones_particles[i].startswith(" Nparticles"):
            particles = np.array(
                zones_particles[i + 2].split(" ")[:-1], dtype=int
            )
            particles_in_zones[f"{particles[0]}"] = particles  # <--
    return particles_in_zones


def create_zobov_voids_properties_and_particles(
    *, jozov_text_file_output_path, particle_in_voids
):
    """
    Create VoidProperties objects from ZOBOV output and particle data.

    Parameters
    ----------
    jozov_text_file_output_path : str
        Path to the JOZOV text file output.
    particle_in_voids : dict
        Dictionary of particles in voids from _get_particles_in_void.

    Returns
    -------
    tuple
        Tuple of (VoidProperties, particles) pairs for voids found by ZOBOV.
    """
    with open(jozov_text_file_output_path, "r") as f:
        voids = f.readlines()

    parameters = [
        "void_number",
        "file_void_number",
        "core_particle",
        "core_dens",
        "zone_vol",
        "zone_number_part",
        "void_number_zones",
        "void_vol",
        "void_number_part",
        "void_dens_contrast",
        "void_prob",
    ]

    properties_and_particles = []
    for void in voids[2:]:
        properties = dict(zip(parameters, void.split()))
        particles = particle_in_voids[str(properties["core_particle"])]
        zobov_void_properties = VoidProperties(**properties)
        properties_and_particles.append((zobov_void_properties, particles))

    return tuple(properties_and_particles)


def process_and_extract_void_properties_and_particles(
    *,
    tinz_output_file_path,
    zinv_output_file_path,
    jozov_text_file_output_path,
):
    """Process ZOBOV output files and extract void properties and particles.

    This function performs the following steps:
    1. Parses tracers in zones raw file.
    2. Creates an array of particles for each void.
    3. Generates VoidProperties objects for each void.

    Parameters
    ----------
    tinz_output_file_path : str
        Path where the parsed tracers vs zones output will be saved.
    zinv_output_file_path : str
        Path where the parsed zones vs voids output will be saved.
    jozov_text_file_output_path : str
        Path to the JOZOV text file output.

    Returns
    -------
    tuple
        Tuple of (VoidProperties, particles) pairs for voids found by ZOBOV.
    """
    # Create dictionary of array of particles
    p_in_v = _get_particles_in_void(
        txt_path=jozov_text_file_output_path,
        tracers_in_zones_path=tinz_output_file_path,
        zones_in_void_path=zinv_output_file_path,
    )

    # Get list of ZobovVoids objects => [(ZVP, PTS), (ZVP, PTS)]
    void_properties_and_particles = (
        create_zobov_voids_properties_and_particles(
            jozov_text_file_output_path=jozov_text_file_output_path,
            particle_in_voids=p_in_v,
        )
    )

    return void_properties_and_particles


def get_zones_in_void(zones_in_void_file_path):
    """
    Gets zones belonging to voids.

    Read the output file containing zones in each void and returns an array
    maping zones to the void they belong.

    Parameters
    ----------
    zones_in_void_file_path : str
        Path to the output file containing zones in each void

    Returns
    -------
    list of numpy.ndarray
        A list of numpy arrays where the first element of each array is an 
        index. The following elements are the zones inside the void, with 
        the void index being the same as the first element of the array.
    """
    with open(zones_in_void_file_path, "r") as f:
        zones = f.readlines()
    zones_in_void = [np.array(zone.split(), dtype=int) for zone in zones[2:]]
    return zones_in_void


def _get_file_void_and_core_particle(*, txt_path):
    """
    Extract void and core particle data from a specified text file.

    Returns the FileVoid# and CoreParticle of each void.

    Parameters
    ----------
    txt_path : str
        Path to the text file containing void and core particle data.

    Returns
    -------
    numpy.ndarray
        An array where each element contains the FileVoid# and 
        CoreParticle associated with each void.
    """
    fv_cp = []
    with open(txt_path, "r") as f:
        data = f.readlines()
        for line in data[2:]:
            line_split = line.split()
            fv_cp.append(np.array([line_split[1], line_split[2]], dtype=int))
    fv_cp = np.array(fv_cp)
    return fv_cp


def _get_particles_in_void(
    *, txt_path, tracers_in_zones_path, zones_in_void_path
):
    """
    Create an array of tracers located inside each void.

    Creates an array of tracers inside each void.

    Parameters
    ----------
    txt_path : str
        Path to the ASCII text file resulting from JOZOV which contains the
        properties of the void.
    tracers_in_zones_path : str
        Path to the file with particles in zones data.
    zones_in_void_path : str
        Path to the file that contains the zones in void data.

    Returns
    -------
    dict
        A dictionary where each key is the core particle of a void, and 
        each value is a numpy array of integers referencing particles 
        in the input box. The first particle in each array is the core 
        particle of the corresponding void (see VoidProperties).
        
    Notes
    -----
    This function may be slow due to its concatenation operations.
    """
    # Get the tracers in zones from the parsed file
    tracers_in_zones = _get_particles_in_zones(
        particles_in_zones_path=tracers_in_zones_path
    )
    # Get the zones in each void from the parsed file
    zones_in_void = get_zones_in_void(
        zones_in_void_file_path=zones_in_void_path
    )
    # Get the columns FileVoid# y CoreParticle
    fv_cp = _get_file_void_and_core_particle(txt_path=txt_path)

    # This part is SLOW
    particles_in_void = {}

    for zones in zones_in_void:

        array = np.array([], dtype=int)
        for zone in zones[1:]:
            index = np.where(fv_cp[:, 0] == zone)[0][0]
            core_particle = fv_cp[index][1]
            array = np.concatenate(
                (array, tracers_in_zones[str(core_particle)]), dtype=int
            )
        particles_in_void[str(array[0])] = array
    return particles_in_void


def get_void_xyz_centers(*, box, txt_path):
    """
    Retrieve xyz coordinates of the CoreParticle for each void.

    Finds xyz coordinates of the CoreParticle of each void.

    Parameters
    ----------
    box : object
        Box object with the tracers properties.
    txt_path : str
        Path to the JOZOV text file containing void properties.

    Returns
    -------
    np.ndarray
        An array of [x, y, z] elements for each CoreParticle.
    """
    xyz = np.column_stack((box.x.value, box.y.value, box.z.value))
    fv_cp = _get_file_void_and_core_particle(txt_path=txt_path)
    fv_cp = fv_cp.T
    return xyz[fv_cp[1]]
