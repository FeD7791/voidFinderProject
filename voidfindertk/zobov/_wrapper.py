import sh

def run_zobov_void_finder(
    path_src, 
    path_input_file,
    buffer_size, 
    box_size, 
    number_of_divisions, 
    executable_name,
    output_name_particles_in_zones,
    output_name_zones_in_void,
    output_name_text_file,
    density_threshold):
    
    # Runing Zobov
    zobov = sh.Command("vozinit",search_paths=[path_src]) #path_src --> pathlib.Path
    zobov(f"{path_input_file}",
          f"{buffer_size}",
          f"{box_size}",
          f"{number_of_divisions}",
          f"{executable_name}")
    # # Run src<executable_name>
    # sh.Command(path_src/f"scr{executable_name}")() 
    
    # # Run jozov
    # jozov = sh.Command(path_src/"jozov")
    # jozov(f"adj{executable_name}.dat",
    #       f"vol{executable_name}",
    #       f"{output_name_particles_in_zones}.dat",
    #       f"{output_name_zones_in_void}.dat",
    #       f"{output_name_text_file}.dat"
    #       f"{density_threshold}")      