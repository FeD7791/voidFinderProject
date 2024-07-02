"""
This module contains several modules parse output files
obtained through the implemetation of the ZOBOV finder algorythm
"""

import ctypes

from astropy import units as u

import numpy as np

import uttr


@uttr.s(repr=False)
class VoidProperties:
    """
    Class that represents the properties of the output ascii
    file of ZOBOV. This Class contains the properties of the
    voids finded with ZOBOV
        Void#:
                What rank the void has, in decreasing order of
                VoidDensContrast.
        FileVoid#:
                The number of this void (starting with 0)
                in the first two files.
        CoreParticle:
                The particle number (starting with 0) of the void's (and
                zone's) core particle (i.e. if CoreParticle=2, it
                would be the third particle in vol...dat and adj...dat).
                Currently, because jozov does not load in particle
                positions (saving memory), looking up the coordinates of
                this particle in the original position file is currently
                the only way of finding the x,y,z position of the core
                of the void. With some prodding, the author might change
                this.
        CoreDens:
                The density, in units of the mean, of the void's
                core particle.
        ZoneVol:
                The volume of the central zone of the void,
                in units of the volume occupied by a mean-density particle.
        Zone#Part:
                The number of particles in the central zone of the void.
        Void#Zones:
                The number of zones in the void.
        VoidVol:
                The volume of the void, in units of the volume occupied by
                a mean-density particle.
        Void#Part:
                The number of particles in the void.
        VoidDensContrast:
                The density contrast of the void, i.e. the ratio between
                the critical density at which water in that zone would flow
                into a deeper zone to the minimum density.
        VoidProb:
                The probability that that DensContrast would arise
                from Poisson noise, using eq. 1 of the ZOBOV paper.
                This probability is based on a fit to the probability
                distribution of DensContrasts from a Poisson particle
                distribution.
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
        """Representation method.

        Returns
        -------
            str
                Name plus number voids found
        """
        cls_name = type(self).__name__
        return f"<{cls_name} void_number={self.void_number}>"


def parse_tracers_in_zones_output(
    *, executable_path, input_file_path, output_file_path
):
    """
    Parse the output raw file from Zobov's output of tracers inside zones
    resulting from the execution of jozov executable.

    Parameters
    ----------
    executable_path : str
        The path to the executable file.
    input_file_path : str
        The path to the input raw file.
    output_file_path : str
        The path to the output ascii file.

    Returns
    -------
    None
        This function does not return any value. It parses the input raw file
        and generates an output raw file.

    Notes
    -----
    This function uses ctypes to interact with the C library generated from
    Zobov's output processing code. The output Ascii file has the following
    format:
    Line 2-6: Total number of particles and total number of voids
    Lines 9 -> : List of the particles with the following format:
        ------------------------
        Nparticles : (Number of particles)
        particulas
        [List of particles]
    First memeber of [List of particles] is the core particle (See
    VoidProperties)
    """
    # Get library
    clibrary = ctypes.CDLL(str(executable_path), mode=ctypes.RTLD_GLOBAL)

    # Input argtypes
    clibrary.get_tracers_in_zones.argtypes = [ctypes.c_char_p, ctypes.c_char_p]

    # Call function
    clibrary.get_tracers_in_zones(
        str(input_file_path).encode(), str(output_file_path).encode()
    )


def parse_zobov(*, jozov_text_file_output_path):
    """
    This method will parse the output of ZOBOV's ascii file into an object:
    VoidProperties (See VoidProperties Class within this module) with the void
    properties.
    Parameters
    ----------
        jozov_text_file_output_path: str
            The path to the file to be parsed. The name is referenced to the
            txt output file that results from the ZOBOV's jozov step.
            This file contains:
            1. Line : Number of particles and voids found by the method.
            2. Line : Name of the attributes asigned to the voids. Attributes
            are:
            -Void# : What rank the void has, in decreasing order of
                    VoidDensContrast.
            -FileVoid# : The number of this void (starting with 0) in the
                        first two files.
            -CoreParticle : The particle number (starting with 0) of the
                            void's (and zone's) core particle (i.e. if
                            CoreParticle=2, it would be the third particle in
                            vol...dat and adj...dat). Currently, because jozov
                            does not load in particle positions (saving memory
                            ), looking up the coordinates of this particle in
                            the original position file is currently the only
                            way of finding the x,y,z position of the core of
                            the void. With some prodding, the author might
                            change this.
            -CoreDens : The density, in units of the mean, of the void's core
                        particle.
            -ZoneVol : The volume of the central zone of the void, in units of
                        the volume occupied by a mean-density particle.
            -Zone#Part : The number of particles in the central zone of the
                        void.
            -Void#Zones :  The number of zones in the void.
            -VoidVol : The volume of the void, in units of the volume occupied
                        by a mean-density particle.
            -Void#Part : The number of particles in the void.
            -VoidDensContrast : The density contrast of the void, i.e. the
                                ratio between the critical density at which
                                water in that zone would flow into a deeper
                                zone to the minimum density.
            -VoidProb : The probability that that DensContrast would arise
                        from Poisson noise, using eq. 1 of the ZOBOV paper.
                        This probability is based on a fit to the probability
                        distribution of DensContrasts from a Poisson particle
                        distribution.
    Returns
    -------
        zobov_voids : list
            List of objects VoidProperties (See VoidProperties whithin
            _methods)
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
        z_void = VoidProperties(**dict(zip(parameters, void.split())))
        zobov_voids.append(z_void)
    return zobov_voids


def get_particles_in_voids(*, particles_in_zones_path):
    """
    This method is used to extract the particles inside each void from the
    parsed file got from the membership file of zobov (which contains
    the particles in each zone). In the parsed file exists n lists of partices
    where the first item in each list is the core particle of the 'core' zone
    of the void. This particle is the one that appears as the CoreParticle in
    the txt output of ZOBOV with the void properties.

    The parsed file is obtained from the ZOBOV's membership file using the
    parse_tracers_in_zones_output private method in this module.

    Parameters
    ----------
        particles_in_zones_path : str
            Path to the parsed file

    Returns
    -------
        dict
            A dictionary where each key is the CoreParticle (see VoidProperties
            in this module) and the values are the list of the particles
            inside each void.
    Notes
    -----
        Please note that each element is not directly related to a void but
        rather is linked indirectly using the CoreParticle as identifier that
        can relate to the void using the txt property file of ZOBOV.
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


def data_post_processing(
    *,
    executable_path,
    input_file_path,
    output_file_path,
    jozov_text_file_output_path,
):
    """
    Process data from the resulting files of a particular run using the ZOBOV
    void finder algorithm.

    Parameters
    ----------
        executable_path: str
            location of the binary file used to parse the code. The default
            binary name is : "tracers_in_zones.so" and is located into
            zobov/src folder
        input_file_path: str
            location of the raw file to be parsed.
        output_file_path: str
            location where the parsed ascii file is going to be located.
        jozov_text_file_output_path: str
            location of the jozov txt output file.
    Returns
    -------
    tuple
        A tuple containing two elements:
        - particle_by_voids : tuple
            Tuple of arrays, where each array contains particles within a void.
        - zobov_voids : list
            List of ZobovVoids objects representing voids parsed from JOZOV
            output.


    """
    # Process 1: Parse tracers in zones raw file in the work directory
    parse_tracers_in_zones_output(
        executable_path=executable_path,
        input_file_path=input_file_path,
        output_file_path=output_file_path,
    )
    # Process 2: Create darray of array of particles
    p_in_v = get_particles_in_voids(particles_in_zones_path=output_file_path)
    # Get list of ZobovVoids objects
    zobov_voids = parse_zobov(
        jozov_text_file_output_path=jozov_text_file_output_path
    )
    # Create Voids as array of particles in void
    voids = []
    for void in zobov_voids:
        voids.append(p_in_v[str(void.core_particle)])
    particle_by_voids = tuple(voids)
    return particle_by_voids, zobov_voids
