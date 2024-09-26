#-----------------------------------------------------------------------------
# COMPILERS

# Select your system: `Linux`, `macOS`, `Host`
# For macOS users, we recommend to install the `gcc` and `hdf5` libraries with Homebrew (https://brew.sh/).
# Check next section for pre-defined Hosts.

SYSTEM = Linux

ifeq ($(SYSTEM),Linux)
    CXX = g++
    MPICXX = mpic++

    #HDF5_HOME =

    CXXFLAGS = -Wall -Wextra -O3 -march=native -fopenmp -std=c++17
    #CXXFLAGS := $(CXXFLAGS) -I$(HDF5_HOME)/include/ -L$(HDF5_HOME)/lib/ -Wl,-rpath,$(HDF5_HOME)/lib
endif

ifeq ($(SYSTEM),macOS)
    CXX = g++-14 # check your version of gcc
    export OMPI_CXX = g++-14
    MPICXX = mpic++

    HDF5_HOME = /opt/homebrew/opt/hdf5
    
    CXXFLAGS = -Wall -Wextra -O3 -march=native -fopenmp -std=c++17
    CXXFLAGS := $(CXXFLAGS) -I$(HDF5_HOME)/include/ -L$(HDF5_HOME)/lib/ -Wl,-rpath,$(HDF5_HOME)/lib
endif

#-----------------------------------------------------------------------------
# PRE-DEFINED HOSTS

# You can add here your pre-defined host.

ifeq ($(SYSTEM),Host)
    HOST=$(shell hostname)
endif

#-- Debugging options:

DEBUG_OPS=-O0 -g -ggdb -fno-inline-small-functions #-gdwarf 

#-- Dante's computer

ifeq ($(HOST),wepwawet)
    CXX = dpcpp
endif

#-- Carlos' Mac

ifeq ($(HOST),op20lm)
    CXX = g++-14
    export OMPI_CXX = g++-14
    MPICXX = mpic++

    HDF5_HOME = /opt/homebrew/opt/hdf5
    
    CXXFLAGS = -Wall -Wextra -O3 -march=native -fopenmp -std=c++17
    CXXFLAGS := $(CXXFLAGS) -I$(HDF5_HOME)/include/ -L$(HDF5_HOME)/lib/ -Wl,-rpath,$(HDF5_HOME)/lib
endif

#-- Clemente cluster (IATE):

ifeq ($(HOST),clemente)
    CXX = g++ # $(DEBUG_OPS)
    MPICXX = mpic++
    
    CXXFLAGS = -Wall -Wextra -O3 -march=native -fopenmp -std=c++17
endif

#-- Raven cluster (MPCDF)

ifeq ($(HOST),raven01)
    CXX = g++ # $(DEBUG_OPS)
    MPICXX = mpic++

    CXXFLAGS = -Wall -Wextra -O3 -march=native -fopenmp -std=c++17
    CXXFLAGS :=$(CXXFLAGS) -I$(HDF5_HOME)/include/ -L$(HDF5_HOME)/lib/ -Wl,-rpath,$(HDF5_HOME)/lib
endif
ifeq ($(HOST),raven02)
    CXX = g++ # $(DEBUG_OPS)
    MPICXX = mpic++

    CXXFLAGS = -Wall -Wextra -O3 -march=native -fopenmp -std=c++17
    CXXFLAGS :=$(CXXFLAGS) -I$(HDF5_HOME)/include/ -L$(HDF5_HOME)/lib/ -Wl,-rpath,$(HDF5_HOME)/lib
endif

#-- 

ifeq ($(HOST),login7a.pri.cosma7.alces.network)
    CXX = g++ # $(DEBUG_OPS)
    MPICXX = mpic++

    CXXFLAGS =-fopenmp  #-I/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/include
endif

#-----------------------------------------------------------------------------
# BUILDING PROCESS

CXXFLAGS +=
LIBS =  -lhdf5
VAR_DIR = configuration_lib
IO_DIR = io_lib
VAR_SOURCE := variables_manager.cpp variables_loader.cpp string_tools.cpp \
			  converter.cpp variables_tools.cpp variable_control.cpp

IO_SOURCE := io.cpp io_tracers_ascii.cpp io_tracers_gadget.cpp \
			  io_tracers_hdf5.cpp io_tracers_stream.cpp

LIB_POP_SOURCE := t.cpp allvars.cpp arvo_variables_tools.cpp arvo_tools.cpp \
	              arvo_functions.cpp objects.cpp grid.cpp

SVF_SOURCE := allvars.cpp t.cpp arvo_variables_tools.cpp arvo_tools.cpp \
	          arvo_functions.cpp objects.cpp grid.cpp finder.cpp svf.cpp

OBJECTS0 := svf.o
OBJECTS1 := popcorn.o
OBJECTS2 := clean_duplicates.o
OBJECTS3 := compute_intersecs.o

OBJECTS_LIB_VAR = $(VAR_SOURCE:.cpp=.o)
OBJECTS_LIB_IO = $(IO_SOURCE:.cpp=.o)
OBJECTS_LIB_DIR = $(addprefix $(VAR_DIR)/, $(OBJECTS_LIB_VAR)) \
                  $(addprefix $(IO_DIR)/, $(OBJECTS_LIB_IO))

OBJECTS_LIB_POP := $(LIB_POP_SOURCE:.cpp=.o)
OBJECTS_LIB_SVF := $(SVF_SOURCE:.cpp=.o)

OBJS = allvars.o t.o arvo_variables_tools.o arvo_tools.o arvo_functions.o \
       objects.o grid.o finder.o svf.o\

INCL = allvars.h Makefile
VAR_TEST = $(VAR_DIR)/unittest
OBJECTS_SVF := $(SVF_SOURCE:.cpp=.o)

TARGETS=svf popcorn_app popcorn clean_duplicates compute_intersecs

all: $(TARGETS)

$(OBJS): $(INCL)

.PHONY: all

clean:
	@echo cleaning bin/
	rm -f $(OBJECTS_LIB_DIR) $(OBJECTS_LIB_POP)  $(VAR_TEST)
	rm -f *.o
	rm -f $(TARGETS)
.PHONY: clean

popcorn_app: popcorn_app.o $(OBJECTS0) $(OBJECTS1) $(OBJECTS2) $(OBJECTS3) $(OBJS) $(OBJECTS_LIB_DIR) $(OBJECTS_LIB_POP)
	$(CXX) $(CXXFLAGS) $(LINKFLAGS) $^ $(LIBS) -o $@

svf: $(OBJECTS_LIB_DIR) $(OBJECTS1) comm.o $(OBJS) $(OBJECTS0:.o=_run.o)
	$(MPICXX) $(CXXFLAGS) $(LINKFLAGS) $^ $(LIBS) -o $@

# popcorn: $(OBJECTS1:.o=_run.o)  $(OBJECTS1)  $(OBJECTS_LIB_DIR)  $(OBJECTS_LIB_POP)
#	$(CXX) $(CXXFLAGS) $(LINKFLAGS) $^ $(LIBS) -o $@

popcorn: $(OBJECTS1:.o=_run.o) $(OBJS) comm.o $(OBJECTS1)  $(OBJECTS_LIB_DIR)  $(OBJECTS_LIB_POP)
	$(MPICXX) $(CXXFLAGS) $(LINKFLAGS) $^ $(LIBS) -o $@

clean_duplicates: $(OBJECTS2:.o=_run.o) $(OBJECTS2) $(OBJECTS_LIB_DIR) $(OBJECTS_LIB_POP)
	$(CXX) $(CXXFLAGS) $(LINKFLAGS) $^ $(LIBS) -o $@

compute_intersecs: $(OBJECTS3:.o=_run.o) $(OBJECTS3) $(OBJECTS_LIB_DIR) $(OBJECTS_LIB_POP)
	$(CXX) $(CXXFLAGS) $(LINKFLAGS) $^ $(LIBS) -o $@

popcorn_run.o: popcorn_run.cpp
	$(MPICXX) $(CXXFLAGS) $(LIBS) -c $< -o $@

svf_run.o: svf_run.cpp
	$(MPICXX) $(CXXFLAGS) $(LIBS) -c $< -o $@

comm.o: comm.cpp
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

#$(OBJECTS_LIB_VAR): $(addprefix $(VAR_DIR)/, $(VAR_SOURCE))
#%.o: %.cpp $(INCLUDES)
#	$(CXX) $(CXXFLAGS) -c $< -o $@
