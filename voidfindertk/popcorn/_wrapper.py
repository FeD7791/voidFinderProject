import sh
import configparser
from ..utils import chdir

# Reference:
# https://gitlab.com/dante.paz/popcorn_void_finder#43-popcorn-void-finder


def popcorn_void_finder(*,mpi_flags,bin_path, conf_file_path, work_dir_path):
    popcorn = sh.Command("popcorn",search_paths=[bin_path])
    params = "config=" + str(conf_file_path)
    # Command will be executed from work_dir_path path.
    with chdir(work_dir_path):
        output = popcorn(params)
    return output

def compute_intersects(*,bin_path, conf_file_path, work_dir_path):
    compute_intersecs = sh.Command("compute_intersecs",search_paths=[bin_path])
    params = "config=" + str(conf_file_path)
    # Command will be executed from work_dir_path path.
    with chdir(work_dir_path):
        output = compute_intersecs(params)
    return output

def clean_duplicates(*,bin_path, conf_file_path, work_dir_path):
    clean_duplicates = sh.Command("clean_duplicates",search_paths=[bin_path])
    params = "config=" + str(conf_file_path)
    # Command will be executed from work_dir_path path.
    with chdir(work_dir_path):
        output = clean_duplicates(params)
    return output



def read_and_modify_config(*,config_file_path, section, parameter, new_value):
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    config.optionxform = str
    # Read the configuration file
    config.read(config_file_path)

    # Check if the section exists
    if not config.has_section(section):
        print(f"Section '{section}' not found in the configuration file.")
        return

    # Check if the parameter exists
    if not config.has_option(section, parameter):
        print(f"Parameter '{parameter}' not found in section '{section}'.")
        return

    # Modify the parameter value
    config.set(section, parameter, new_value)

    # Save the changes back to the configuration file
    with open(config_file_path, 'w') as configfile:
        config.write(configfile)