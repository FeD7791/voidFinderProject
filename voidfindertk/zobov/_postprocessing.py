"""This module contains functions to parse output files from the ZOBOV \
void finder.

"""

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


def get_particles_in_voids(*, particles_in_zones_path):
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
            particles_in_zones[f"{particles[0]}"] = particles
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
        Dictionary of particles in voids from get_particles_in_voids.

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
    executable_path,
    input_file_path,
    output_file_path,
    jozov_text_file_output_path,
):
    """Process ZOBOV output files and extract void properties and particles.

    This function performs the following steps:
    1. Parses tracers in zones raw file.
    2. Creates an array of particles for each void.
    3. Generates VoidProperties objects for each void.

    Parameters
    ----------
    executable_path : str
        Path to the C library executable for parsing tracers.
    input_file_path : str
        Path to the input file with tracers in zones data.
    output_file_path : str
        Path where the parsed tracers output will be saved.
    jozov_text_file_output_path : str
        Path to the JOZOV text file output.

    Returns
    -------
    tuple
        Tuple of (VoidProperties, particles) pairs for voids found by ZOBOV.
    """
    # Process 1: Parse tracers in zones raw file in the work directory
    parse_tracers_in_zones_output(
        executable_path=executable_path,
        input_file_path=input_file_path,
        output_file_path=output_file_path,
    )

    # Process 2: Create dictionary of array of particles
    p_in_v = get_particles_in_voids(particles_in_zones_path=output_file_path)

    # Get list of ZobovVoids objects => [(ZVP, PTS), (ZVP, PTS)]
    void_properties_and_particles = (
        create_zobov_voids_properties_and_particles(
            jozov_text_file_output_path=jozov_text_file_output_path,
            particle_in_voids=p_in_v,
        )
    )

    return void_properties_and_particles
