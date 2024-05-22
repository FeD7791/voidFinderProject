#Compiler
CC = g++


##Paths Voro++
### -L[/path/to/shared-libraries]
### -I[/path/to/header-files]

VORO_L = ./voro++-0.4.6/src
VORO_I = ./voro++-0.4.6/src

##Paths GSL
GSL_L = /opt/spack/dev/opt/spack/linux-centos7-x86_64_v3/gcc-12.1.0/gsl-2.7.1-ttp5l7ez636aw4kkpbfzfufdu2546jgf/lib
GSL_I = /opt/spack/dev/opt/spack/linux-centos7-x86_64_v3/gcc-12.1.0/gsl-2.7.1-ttp5l7ez636aw4kkpbfzfufdu2546jgf/include

#OPEN MP flag
OMP = -fopenmp

#GSL Flags
GSL_F = -lgsl -lgslcblas

#VORO++ Flags
VORO_F = -lvoro++

#WARNING ALL
WARNINGS = -Wall

##CCFILES
CCF = pipeline_vf.cc tools.cc cosmology.cc io.cc voronoi.cc grid.cc finder.cc qsort.cc velocity.cc profiles.cc allvars.cc
##CCOBJETS
CCO = pipeline_vf.o tools.o cosmology.o io.o voronoi.o grid.o finder.o qsort.o velocity.o profiles.o allvars.o


#Compilation step 1
.PHONY: install
install:
	$(CC) -c $(WARNINGS) $(OMP) -fPIC $(CCF) -I $(VORO_I) -I $(GSL_I) -I $(GSL_I)/gsl
	$(CC) -shared -o vf_lib.so $(OMP) $(CCO) -L $(GSL_L) -L $(VORO_L) -Wl,--whole-archive  $(GSL_F) $(VORO_F) -lm -Wl,--no-whole-archive

#LD_LIBRARY_PATH (if necesary)
.PHONY: export_lib
export_lib:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'$(GSL_L)'

#Istall Voro++
.PHONY: install_voro++
install_voro++:
	wget https://math.lbl.gov/voro++/download/dir/voro++-0.4.6.tar.gz
	tar -xvf  voro++-0.4.6.tar.gz
	cp config.mk ./voro++-0.4.6
	cd ./voro++-0.4.6 && $(MAKE)

.PHONY: clean
clean:
	rm -f $(CCO) vf_lib.so
