from configparser import ConfigParser
import pandas as pd
import sh
from ..utils import chdir

def config_file_maker(
    *,
    # Parameters
    trsfile,
    filefmt,
    num_file,
    sphfile,
    popfile,
    auxfiles,
    rawpopfile,
    pairsfile,
    boxsize,
    densth,
    minradius,
    maxradius,
    massmin,
    eps,
    # Directory path to place the file
    path
    ):
    '''
    Generates the configuration file with the parameters to run the void
    finder in the desired path.
    Parameters
    ----------
    trsfile(str): 
        Input tracer file.
    filefmt(str):
        File format options: "ASCII" "STREAM" "HDF5" "HDF5_SUBFIND_GROUPS"
        "HDF5_SUBFIND_SUBHALOS" "GADGET1" "GADGET2" "GADGET4_TYPE1"
    num_file(str):
        Number of files for the tracer catalogue, specially important for
        Gadget outputs.
    sphfile(str):
        Output Spherical voids catalogue.
    popfile(str):
        Output Popcorn voids catalogue (after cleaning overlappings).
    auxfiles(str):
        Whether or not use auxilary files:{true , false}
    rawpopfile(str):
        File with Popcorn voids before cleaning overlappings.
    pairsfile(str):
        File with Pairs of touching Popcorn voids.
    boxsize(str):
        Length of the box in the same units of the tracer input file
        (required only for ascii inputs, otherwise ignored).
    densth(str):
        Integrated density threshold for identification.
    minradius(str):
        Minimal radii allowed for a sphere member in input units.
    maxradius(str):
        Maximal radii allowed for a sphere member in input units.
    massmin(str):
        Minimal halo mass allowed (set to 0 if not applicable).
    eps(str):
        Obsolete flag, don't modify.
    path(pathlib.Path):
        Path where the file is going to be generated.
    '''
    config = ConfigParser(allow_no_value=True)
    config.optionxform = str
    config.add_section('INPUT_PARAMS')
    config.set('INPUT_PARAMS', 'TRSFILE', trsfile)
    config.set('INPUT_PARAMS', 'FILEFMT', filefmt)
    config.set('INPUT_PARAMS', 'NUM_FILE', num_file)
    config.set('INPUT_PARAMS', 'SPHFILE', sphfile)
    config.set('INPUT_PARAMS', 'POPFILE', popfile)
    config.set('INPUT_PARAMS', 'AUXFILES', auxfiles)
    config.set('INPUT_PARAMS', 'RAWPOPFILE', rawpopfile)
    config.set('INPUT_PARAMS', 'PAIRSFILE', pairsfile)
    config.set('INPUT_PARAMS', 'BOXSIZE', boxsize)
    config.set('INPUT_PARAMS', 'DENSTH', densth)
    config.set('INPUT_PARAMS', 'MINRADIUS',minradius)
    config.set('INPUT_PARAMS', 'MAXRADIUS', maxradius)
    config.set('INPUT_PARAMS', 'MASSMIN', massmin)
    config.set('INPUT_PARAMS', 'EPS', eps)

    with open(path, 'w') as configfile:
        config.write(configfile)

def popcorn_svf_input_data_builder(*,box,file_path):
    '''
    Generates input file from box in the desired file path.
    Parameters
    ----------
    box(Box object):
        Box object with the tracers and their properties.
    file_path(pathlib.Path):
        File path to place the file.
    '''
    df = pd.DataFrame(box.__dict__)
    df.drop(labels=['_len'], axis=1,inplace=True)
    df = df[["m","x","y","z","vx","vy","vz"]] #Popcorn input file format
    df.to_csv(
        file_path,
        sep=' ',
        index=False,
        header=False,
        float_format='%.2f'
        )

def spherical_popcorn_void_finder(
    *,
    mpi_flags,
    bin_path,
    conf_file_path,
    work_dir_path
    ):
    '''
    Runs the Popcorn-Spherical Void Finder using mpi.
    Parameters
    ----------
    mpi_flags(str):
        mpi flags.
    bin_path(pathlib.Path):
        Path to the folder where the binary file that runs the void finder is.
    conf_file_path(pathlib.Path):
        Path to the configuration file used to run de void finder.
    work_dir_path : str
        Path to the working directory.
    Returns
    -------
    output:
        sh output of the run.
    '''
    # Reference to mpi
    #https://docs.oracle.com/cd/E19356-01/820-3176-10/ExecutingPrograms.html
    svf_mpi = sh.Command('svf',search_paths=[bin_path])
    params = 'config='+str(conf_file_path / "vars.conf")
    # Command will be executed from work_dir_path path.
    with chdir(work_dir_path):
        output = svf_mpi(params)
    return output
    #subprocess.run(["mpirun", "-np", "1","--bind-to","none", "./svf", "./configuration/vars.conf"])
