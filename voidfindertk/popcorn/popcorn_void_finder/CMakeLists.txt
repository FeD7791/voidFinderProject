cmake_minimum_required(VERSION 3.12)
project(Popcorn)

# Set compiler options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -march=native -fopenmp")

# Find MPI package
find_package(MPI REQUIRED)

# Include directories for HDF5 (if needed)
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Source files
set(VAR_DIR ./Source/configuration_lib)
set(SRC_DIR ./Source)
set(IO_DIR  ./Source/io_lib)

set(VAR_SOURCE ${VAR_DIR}/variables_manager.cpp
               ${VAR_DIR}/variables_loader.cpp
               ${VAR_DIR}/string_tools.cpp
               ${VAR_DIR}/converter.cpp
               ${VAR_DIR}/variables_tools.cpp
               ${VAR_DIR}/variable_control.cpp)

set(IO_SOURCE ${IO_DIR}/io.cpp
              ${IO_DIR}/io_tracers_ascii.cpp
              ${IO_DIR}/io_tracers_gadget.cpp
              ${IO_DIR}/io_tracers_hdf5.cpp
              ${IO_DIR}/io_tracers_stream.cpp)


set(POP_SOURCE ${SRC_DIR}/allvars.cpp
               ${SRC_DIR}/t.cpp
               ${SRC_DIR}/arvo_variables_tools.cpp
               ${SRC_DIR}/arvo_tools.cpp
               ${SRC_DIR}/arvo_functions.cpp
               ${SRC_DIR}/objects.cpp
               ${SRC_DIR}/grid.cpp
               ${SRC_DIR}/comm.cpp
               ${SRC_DIR}/finder.cpp
               ${SRC_DIR}/popcorn.cpp
               ${SRC_DIR}/compute_intersecs.cpp
               ${SRC_DIR}/clean_duplicates.cpp
               ${SRC_DIR}/svf.cpp)

# Create library targets
add_library(var_lib OBJECT ${VAR_SOURCE})
add_library(io_lib OBJECT ${IO_SOURCE})

# Executable targets
add_executable(svf ${SRC_DIR}/svf_run.cpp $<TARGET_OBJECTS:var_lib> $<TARGET_OBJECTS:io_lib> ${POP_SOURCE})
add_executable(popcorn ${SRC_DIR}/popcorn_run.cpp $<TARGET_OBJECTS:var_lib> $<TARGET_OBJECTS:io_lib> ${POP_SOURCE})
add_executable(popcorn_app ${SRC_DIR}/popcorn_app.cpp $<TARGET_OBJECTS:var_lib> $<TARGET_OBJECTS:io_lib> ${POP_SOURCE})
add_executable(clean_duplicates ${SRC_DIR}/clean_duplicates_run.cpp $<TARGET_OBJECTS:var_lib> $<TARGET_OBJECTS:io_lib> ${POP_SOURCE})
add_executable(compute_intersecs ${SRC_DIR}/compute_intersecs_run.cpp $<TARGET_OBJECTS:var_lib> $<TARGET_OBJECTS:io_lib> ${POP_SOURCE})

# Link against HDF5 library
target_link_libraries(svf PRIVATE MPI::MPI_CXX ${HDF5_LIBRARIES})
target_link_libraries(popcorn_app PRIVATE MPI::MPI_CXX ${HDF5_LIBRARIES})
target_link_libraries(popcorn PRIVATE MPI::MPI_CXX ${HDF5_LIBRARIES})
target_link_libraries(clean_duplicates PRIVATE MPI::MPI_CXX ${HDF5_LIBRARIES})
target_link_libraries(compute_intersecs PRIVATE MPI::MPI_CXX ${HDF5_LIBRARIES})

# Clean target
add_custom_target(clean-all
        COMMAND ${CMAKE_BUILD_TOOL} clean
        COMMAND rm -f CMakeCache.txt cmake_install.cmake Makefile
        COMMAND rm -rf CMakeFiles/
        COMMAND rm -rf ${VAR_DIR}/*.o ${IO_DIR}/*.o
        COMMAND rm -f *.o popcorn_app svf popcorn clean_duplicates compute_intersecs
        )

# Set output directories
set_target_properties(svf popcorn_app popcorn clean_duplicates compute_intersecs
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        )
