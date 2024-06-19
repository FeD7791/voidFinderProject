import ctypes

from astropy import units as u

import numpy as np

import pandas as pd

import uttr


@uttr.s(repr=False)
class _ZobovVoids:

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

    Void_number = uttr.ib(converter=np.array)
    File_void_number = uttr.ib(converter=np.array)
    CoreParticle = uttr.ib(converter=np.array)
    CoreDens = uttr.ib(converter=np.array)
    ZoneVol = uttr.ib(converter=np.array)
    Zone_number_part = uttr.ib(converter=np.array)
    Void_number_Zones = uttr.ib(converter=np.array)
    VoidVol = uttr.ib(converter=np.array, unit=u.Mpc**3)
    Void_number_Part = uttr.ib(converter=np.array)
    VoidDensContrast = uttr.ib(converter=np.array)
    VoidProb = uttr.ib(converter=np.array)
    _void_len = uttr.ib(init=False)
    _tracers_in_void = uttr.ib(
        init=False, default=None
    )  # Provide tracers in void

    def __attrs_post_init__(self):
        """Post init method.

        Checks that the lenght of the inputs are the same
        """
        lengths = set()
        for e in (
            # Rank of the void, in decreasing order of VoidDensContrast.
            self.Void_number,
            # Number of the void, in previous files
            self.File_void_number,
            # Index --Related to Box-- of the core particle of the void
            self.CoreParticle,
            # Density of the core particle
            self.CoreDens,
            # Volume of the central zone of the void
            self.ZoneVol,
            # Number of particles in the central zone of the void
            self.Zone_number_part,
            # Number of zones in the void
            self.Void_number_Zones,
            # Volume of the void
            self.VoidVol,
            # Number of particles in the void
            self.Void_number_Part,
            # Contrast density of the void
            self.VoidDensContrast,
            # The Poisson noise probability
            self.VoidProb,
        ):
            lengths.add(len(e))

        if len(lengths) != 1:
            raise ValueError("Arrays should be of the same size")

        super().__setattr__("_void_len", lengths.pop())

    def __len__(self):
        """Length method.

        Returns
        -------
            int
                the number of elements in SphericalVoids
        """
        return self._void_len

    def __repr__(self):
        """Representation method.

        Returns
        -------
            str
                Name plus number of points in the box
        """
        cls_name = type(self).__name__
        length = len(self)
        return f"<{cls_name} size={length}>"


def parse_zones_in_voids_output(
    *, executable_path, input_file_path, output_file_path
):
    """
    Parse the output raw file from Zobov's output of zones inside voids
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
    # path = _Paths.ZOBOV / "out_zones_in_void.dat"
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


def read_zobov_output(filename_path):
    """
    This method will parse the output of ZOBOV's ascii file
    into an object: ZobovVoids with the void properties
    """
    output = pd.read_csv(filename_path, sep="\s+", skiprows=2)
    output.columns = [
        "Void_number",
        "File_void_number",
        "CoreParticle",
        "CoreDens",
        "ZoneVol",
        "Zone_number_part",
        "Void_number_Zones",
        "VoidVol",
        "Void_number_Part",
        "VoidDensContrast",
        "VoidProb",
    ]
    zobov = _ZobovVoids(**output.to_dict(orient="list"))
    return zobov
