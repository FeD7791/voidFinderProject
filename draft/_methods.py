"""
This module contains several modules parse output files
obtained through the implemetation of the ZOBOV finder algorythm
"""

import ctypes

from astropy import units as u

import numpy as np

import pandas as pd

import uttr


@uttr.s(repr=False)
class ZobovVoids:
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
    _void_len = uttr.ib(init=False)
    tracers_in_void = uttr.ib(
        init=False, default=None
    )  # Provide tracers in void

    def __repr__(self):
        """Representation method.

        Returns
        -------
            str
                Name plus number voids found
        """
        cls_name = type(self).__name__
        return f"<{cls_name} void_number={self.void_number}>"


# def _parse_zones_in_voids_output(
#     *, executable_path, input_file_path, output_file_path
# ):
#     """
#     Parse the output raw file from Zobov's output of zones inside voids
#     resulting from the execution of jozov executable.

#     Parameters
#     ----------
#     executable_path : str
#         The path to the executable file.
#     input_file_path : str
#         The path to the input raw file.
#     output_file_path : str
#         The path to the output ascii file.

#     Returns
#     -------
#     None
#         This function does not return any value. It parses the input raw file
#         and generates an output Ascii file.

#     Notes
#     -----
#     This function uses ctypes to interact with the C library generated from
#     Zobov's output processing code.
#     """
#     # path = _Paths.ZOBOV / "out_zones_in_void.dat"
#     # Get library
#     clibrary = ctypes.CDLL(str(executable_path), mode=ctypes.RTLD_GLOBAL)

#     # Input argtypes
#     clibrary.process_files.argtypes = [ctypes.c_char_p, ctypes.c_char_p]

#     # Call function
#     clibrary.process_files(
#         str(input_file_path).encode(), str(output_file_path).encode()
#     )


def _parse_tracers_in_zones_output(
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
    Zobov's output processing code.
    """
    # Get library
    clibrary = ctypes.CDLL(str(executable_path), mode=ctypes.RTLD_GLOBAL)

    # Input argtypes
    clibrary.get_tracers_in_zones.argtypes = [ctypes.c_char_p, ctypes.c_char_p]

    # Call function
    clibrary.get_tracers_in_zones(
        str(input_file_path).encode(), str(output_file_path).encode()
    )


def parse_zobov(filename_path):
    """
    This method will parse the output of ZOBOV's ascii file
    into an object: ZobovVoids with the void properties
    """

    with open(filename_path, "r") as f:
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
        z_void = ZobovVoids(**dict(zip(parameters, void.split())))
        zobov_voids.append(z_void)

    return zobov_voids


# def _get_zones_in_voids(*,zones_in_voids_file_path:str):
#     """
#     Reads zone data from the file containing zones in voids information.
#     This file has the following format:
#         First two rows are : space / number of voids
#         The rest of the filecontest : n_void - list of zones in the void

#     Parameters
#     ----------
#     zones_in_voids_file_path : str
#         Path to the file containing void zones information.

#     Returns
#     -------
#     dict
#         A dictionary where keys are void number and values are NumPy arrays of
#         integers representing zones in that void.
#     """
#     # Read the file and save its contents in a list
#     with open(zones_in_voids_file_path,"r") as f:
#         void_zones = f.readlines()

#     zones_in_voids = {} # Dict to hold the zones in each void
#     for vz in void_zones[2:]:
#         vz = vz.split(" ") #split each element by a space

#         # The void number is vz[0] the zones are contained in vz[1:-1]
#         zones_in_voids[vz[0]] = np.array(vz[1:-1],dtype=int)
#     return zones_in_voids


def _get_particles_in_voids(*, particles_in_zones_path):
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
            A dictionary where each key is the CoreParticle (see ZobovVoids
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

    particles_in_zones = {}
    for i in np.arange(len(zones_particles)):

        if zones_particles[i].startswith(" Nparticles"):

            particles = np.array(
                zones_particles[i + 2].split(" ")[:-1], dtype=int
            )
            particles_in_zones[f"{particles[0]}"] = particles
    return particles_in_zones
