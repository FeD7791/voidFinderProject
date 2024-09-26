# Compiler
CXX = g++

# Compiler flags
FLAGS0 = -std=c++11 -fopenmp
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

# Directory structure (adjust as needed for relative paths)
dirLib = /home/jorgefederico/project/CosmoBolognaLib/CosmoBolognaLib/ # Here your directory to the CBL Library
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.4.0/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2.1/

# Library and include flags
FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

# Object files
OBJ = cleaner.o

# Shared library extension based on OS
ES = so
SYS := $(shell uname -s)

ifeq ($(SYS),Darwin)
    ES = dylib
endif

# Default target
all: cleaner

# Target to build the shared library
cleaner: $(OBJ)
	$(CXX) -shared -fPIC -o libcleaner.$(ES) $(OBJ) $(FLAGS_LIB)

# Clean up build artifacts
clean:
	rm -f *.o libcleaner.$(ES) *~ \#* temp* core*

# Rule to compile .cpp to .o
%.o: %.cpp
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c -o $@ $<

# Rule to compile cleaner.cpp
cleaner.o: cleaner.cpp
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c -fPIC -o cleaner.o cleaner.cpp