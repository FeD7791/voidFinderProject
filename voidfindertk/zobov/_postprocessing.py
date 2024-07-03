"""
This module contains several functions to parse output files
obtained through the implementation of the ZOBOV void finder algorithm.
"""

import ctypes
from astropy import units as u
import numpy as np
import uttr


@uttr.s(repr=False)
class VoidProperties:
    """Class that represents the properties of voids found with ZOBOV.

    This class contains the properties of voids identified by ZOBOV,
    as extracted from the output ASCII file.

    Attributes
    ----------
    void_number : int
        Rank of the void in decreasing order of VoidDensContrast.
    file_void_number : int
        Number of this void (starting with 0) in the first two files.
    core_particle : int
        Particle number (starting with 0) of the void's (and zone's) core
        particle.
    core_dens : float
        Density of the void's core particle, in units of the mean density.
    zone_vol : float
        Volume of the central zone of the void, in units of the volume occupied
        by a mean-density particle.
    zone_number_part : int
        Number of particles in the central zone of the void.
    void_number_zones : int
        Number of zones in the void.
    void_vol : float
        Volume of the void, in units of the volume occupied by a mean-density
        particle.
    void_number_part : int
        Number of particles in the void.
    void_dens_contrast : float
        Density contrast of the void, i.e., the ratio between the critical
        density at which water in that zone would flow into a deeper zone to
        the minimum density.
    void_prob : float
        Probability that the DensContrast would arise from Poisson noise,
        using eq. 1 of the ZOBOV paper.
    particles : numpy.ndarray
        Array of particle IDs in the void.

    Notes
    -----
    The void_prob is based on a fit to the probability distribution of
    DensContrasts from a Poisson particle distribution.

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
    particles = uttr.ib(converter=np.array)

    def __repr__(self):
        """Return a string representation of the VoidProperties object.

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
    """Parse tracers in zones output using a C library.

    Parameters
    ----------
    executable_path : str
        Path to the C library executable.
    input_file_path : str
        Path to the input file containing tracers in zones data.
    output_file_path : str
        Path where the parsed output will be saved.

    Notes
    -----
    This function uses ctypes to call a C function that parses the tracers in
    zones output.

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
    """Extract particles in voids from a parsed file.

    Parameters
    ----------
    particles_in_zones_path : str
        Path to the file containing parsed particles in zones data.

    Returns
    -------
    dict
        Dictionary with zone IDs as keys and numpy arrays of particle IDs as
        values.

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


def create_zobov_voids_properties(
    *, jozov_text_file_output_path, particle_in_voids
):
    """Create VoidProperties objects from ZOBOV output and particle data.

    Parameters
    ----------
    jozov_text_file_output_path : str
        Path to the JOZOV text file output.
    particle_in_voids : dict
        Dictionary containing particles in voids, as returned by
        get_particles_in_voids.

    Returns
    -------
    tuple
        Tuple of VoidProperties objects representing the voids found by ZOBOV.

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

    zobov_voids = []
    for void in voids[2:]:
        properties = dict(zip(parameters, void.split()))
        properties["particles"] = particle_in_voids[
            str(properties["core_particle"])
        ]
        z_void = VoidProperties(**properties)
        zobov_voids.append(z_void)

    return tuple(zobov_voids)


def process_and_extract_void_properties(
    *,
    executable_path,
    input_file_path,
    output_file_path,
    jozov_text_file_output_path,
):
    """Process ZOBOV output files and extract void properties.

    This function performs the following steps:
    1. Parses tracers in zones raw file.
    2. Creates an array of particles for each void.
    3. Generates VoidProperties objects for each void.

    Parameters
    ----------
    executable_path : str
        Path to the C library executable for parsing tracers.
    input_file_path : str
        Path to the input file containing tracers in zones data.
    output_file_path : str
        Path where the parsed tracers output will be saved.
    jozov_text_file_output_path : str
        Path to the JOZOV text file output.

    Returns
    -------
    tuple
        Tuple of VoidProperties objects representing the voids found by ZOBOV.
    """
    # Process 1: Parse tracers in zones raw file in the work directory
    parse_tracers_in_zones_output(
        executable_path=executable_path,
        input_file_path=input_file_path,
        output_file_path=output_file_path,
    )

    # Process 2: Create dictionary of array of particles
    p_in_v = get_particles_in_voids(particles_in_zones_path=output_file_path)

    # Get list of ZobovVoids objects
    void_properties = create_zobov_voids_properties(
        jozov_text_file_output_path=jozov_text_file_output_path,
        particle_in_voids=p_in_v,
    )

    return void_properties
