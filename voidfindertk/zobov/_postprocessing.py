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
    """Class that represents the properties of the output ascii
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
    particles = uttr.ib(converter=np.array)

    def __repr__(self):
        """x.__repr__() <==> repr(x)"""
        cls_name = type(self).__name__
        return f"<{cls_name} void_number={self.void_number}>"


def parse_tracers_in_zones_output(
    *, executable_path, input_file_path, output_file_path
):

    # Get library
    clibrary = ctypes.CDLL(str(executable_path), mode=ctypes.RTLD_GLOBAL)

    # Input argtypes
    clibrary.get_tracers_in_zones.argtypes = [ctypes.c_char_p, ctypes.c_char_p]

    # Call function
    clibrary.get_tracers_in_zones(
        str(input_file_path).encode(), str(output_file_path).encode()
    )


def get_particles_in_voids(*, particles_in_zones_path):

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

    # Process 1: Parse tracers in zones raw file in the work directory
    parse_tracers_in_zones_output(
        executable_path=executable_path,
        input_file_path=input_file_path,
        output_file_path=output_file_path,
    )

    # Process 2: Create darray of array of particles
    p_in_v = get_particles_in_voids(particles_in_zones_path=output_file_path)

    # Get list of ZobovVoids objects
    void_properties = create_zobov_voids_properties(
        jozov_text_file_output_path=jozov_text_file_output_path,
        particle_in_voids=p_in_v,
    )

    # particle_by_voids = tuple(voids)
    return void_properties
