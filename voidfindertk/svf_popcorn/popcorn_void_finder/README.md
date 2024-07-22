# Popcorn Void Finder

We present a new definition of cosmic voids and the code with the algorithm that implements it. Underdense regions are defined as free-form objects, called **popcorn voids**, made from the union of spheres of maximum volume with a given joint <u>integrated (under)density contrast</u>. The method is inspired by the excursion-set theory and consequently no rescaling processing is needed, the removal of overlapping voids and objects with sizes below the shot noise threshold is inherent in the algorithm.

Popcorn voids are an extension of **spherical voids** by adding children spheres in a recursive way. Therefore, the so-called spherical void finder is also included here. Specifically, we provide an enhanced version of the code presented in [Ruiz et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.448.1471R/abstract).

For more details, we invite you to read the following paper: [Paz et al. 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.522.2553P/abstract). If you use this code (or part of it) for a scientific work, we kindly ask you to cite this paper.

<!-- --------------------------------------------------------------------- -->
## 1. Required software

A `c++` compiler and the following libraries:
* `MPI`
* `OpenMP`
* `HDF5`

Optional:
* `GTest`
* `Meson build system`

For **macOS** users, we recommend to install the `gcc` and `hdf5` libraries with **Homebrew** (https://brew.sh/, just `brew install gcc hdf5`).

<!-- --------------------------------------------------------------------- -->
## 2. Installation

We provide two alternative ways to install the code.

### 2.1. Standard installation

* Go to the `Source/` folder.
* Edit the `Makefile` with your preferences. Check the paths to compilers and required libraries. We provide suggestions for **Linux** and **macOS** users, as well as for some pre-defined hosts.
* Compile the code with:
```
make
```

After the building process you should have the following executable binary files:

**Void finders:**
* `svf` :: the spherical mode,
* `popcorn` :: the popcorn mode.

**Auxilary tools:**
* `compute_intersecs` :: to compute all the overlapping popcorn pairs,
* `clean_duplicates` :: to write a clean file removing all small volume overlapping popcorn voids.

### 2.2. Alternative: using Meson

* If you prefer to use the **Meson Build System** (https://mesonbuild.com/), at the `Source/` folder run:
```
meson setup build
```
* After checking dependencies, if everything goes well, in the `build/` folder the project will be ready for compilation:
```
cd build
meson compile
```
### 2.3. Alternative: using CMake

* If you prefer to use the **CMake Build System** (https://cmake.org/), at the principal folder run:
```
 mkdir build
 cd build
 cmake ..
```
* After checking dependencies, if everything goes well, in the `build/` folder the project will be ready for compilation:
```
cmake --build .
```

<!-- --------------------------------------------------------------------- -->
## 3. Preparing your data: input structure

Currently, this version of the Popcorn Void Finder only admits a DM-particle or halo catalogue. We are working on the survey version...

* The allowed formats for particle catalogues are: `GADGET1` (binary), `GADGET2` (binary), `GADGET4_TYPE1` (binary), `HDF5` for Gadget.
* The allowed formats for halo catalogues are: `HDF5_SUBFIND_GROUPS`, `HDF5_SUBFIND_SUBHALOS`.
* Alternatively, you can provide a single `ASCII` file, without header, with the following structure: 

| $n_{\rm part}$ | $x$ | $y$ | $z$ | $v_x$ | $v_y$ | $v_z$ |

**References:**

* $n_{\rm part}$ (`integer`): number of particles of the tracer (e.g., a halo) (set to $0$ if not applicable), 
* $(x,y,z)$ (`float`): cartesian coordinates of the tracer (in arbitrary units, usually in ${\rm Mpc}$ or $h^{-1}{\rm Mpc}$),
* $(v_x,v_y,v_z)$ (`float`): velocity components of the tracer (in ${\rm km}/{\rm s}$). **NOTE:** the velocities are not used for void identification, but they are needed for dynamical analyses that we are planning to incorporate (e.g., calculation of the velocity profile). You can set them to $0$ if not applicable.

<!-- --------------------------------------------------------------------- -->
## 4. Execution

### 4.1. Settings

The file `Source/configuration/vars.conf` contains an example of config file for running a void identification. See the optional parameters and comments there. Here we provide the most basic parameters:
* `BOXSIZE` (`float`): length of the box in the same units of the tracer input file,
* `DENSTH` (`float`): integrated density threshold for identification,
* `MINRADIUS` (`float`): minimal radius allowed for a member sphere in input units,
* `MAXRADIUS` (`float`): maximal radius allowed for a member sphere in input units.

We recommend setting the number of parallel processes (e.g. 8): `export OMP_NUM_THREADS=8`.

If you are using a queuing system (e.g., in a cluster), the `Example/` folder includes a template script for submitting jobs to the **Slurm Workload Manager**.

### 4.2. Spherical void finder

Run the **spherical mode** (in a single node, using openmp):
```
./svf config=<PATH_TO_CONFIG_FILE>/vars.conf
```
A catalogue of spherical voids will be produced.

### 4.3. Popcorn void finder

* Run first the spherical mode, as explained in the previous subsection. Popcorn voids will be obtained starting from these main spheres.
* We recommend adjusting again the `MINRADIUS` parameter in the config file to control the shot-noise level (see our paper for more details).
* Run the **popcorn mode** (in a single or multiple nodes using mpi/openmp):
```
./popcorn config=<PATH_TO_CONFIG_FILE>/vars.conf
```
* Look for superporsitions:
```
./compute_intersecs config=<PATH_TO_CONFIG_FILE>/vars.conf
```
* Clean for overlappings:
```
./clean_duplicates config=<PATH_TO_CONFIG_FILE>/vars.conf
```
A catalogue of popcorn voids will be produced.

<!-- --------------------------------------------------------------------- -->
## 5. Output data structure

### 5.1. Spherical void finder

The output spherical void catalogue is a single `ASCII` file, without header, with the following structure:

| ${\rm ID}$ | $R_{\rm v}$ | $x_{\rm v}$ | $y_{\rm v}$ | $z_{\rm v}$ | $\Delta(R_{\rm v})$ |

**References:**

* ${\rm ID}$ (`integer`): void ID,
* $R_{\rm v}$ (`float`): void radius,
* $(x_{\rm v}, y_{\rm v}, z_{\rm v})$ (`float`): cartesian coordinates of the void centre,
* $\Delta(R_{\rm v})$ (`float`): integrated density contrast at $R_{\rm v}$.

### 5.2. Popcorn void finder

The output popcorn void catalogue is a single `ASCII` file, without header, with the following structure:

* First line: $n_{\rm pop}$.
* Succesive lines: popcorn objets structured as follows:
    * Header line: | ${\rm ID}$ | $n_{\rm mem}$ | ${\rm vol}$ | $n_{\rm part}$ | ${\rm flag}$ |
    * Attributes of member sphere 1: | $x_{\rm v}$ | $y_{\rm v}$ | $z_{\rm v}$ | $R_{\rm v}$ | ${\rm fracvol}$ | ${\rm level}$
    * ... 
    * Attributes of member sphere $n_{\rm mem}$
    * ID of first particle inside the void
    * ...
    * ID of last particle inside the void

**References:**

* $n_{\rm pop}$ (`integer`): number of popcorn voids,
* ${\rm ID}$ (`integer`): popcorn ID,
* $n_{\rm mem}$ (`integer`): number of member spheres,
* ${\rm vol}$ (`float`): popcorn volume,
* $n_{\rm part}$ (`integer`): number of particles inside,
* ${\rm flag}$ (`integer`): for internal control, you can ignore it...,
* Attributes of the member spheres:
    * $x_{\rm v}$ (`float`): x-coordinate of the centre,
    * $y_{\rm v}$ (`float`): y-coordinate,
    * $z_{\rm v}$ (`float`): z-coordinate,
    * $R_{\rm v}$ (`float`): radius,
    * ${\rm fracvol}$ (`float`): volume contribution to the total volume,
    * ${\rm level}$ (`integer`): hierarchy level: 0 for main sphere, 1 for secondary spheres, ...

<!-- --------------------------------------------------------------------- -->
## 6. Running an example

The `Example/` folder contains an example of how to identify voids. We provide a DM-halo catalogue: `halos_ascii.dat`, and a `vars.conf` file suitable for this catalogue. You can download the halo catalogue here: https://iate.oac.uncor.edu/~cmcorrea/halos_ascii.dat.

**Characteristics of the halo catalogue:**

* Format: `ASCII`,
* Origin: A flat-$\Lambda$CDM simulation, $1024^3$ particles, $\Omega_m=0.25$, $h=0.7$, $\sigma_{8/h}=0.9$, snapshot $z=0.51$, in real space,
* Box size: $1000~h^{-1}{\rm Mpc}$,
* Muss cut: $2\times10^{12}~h^{-1}{\rm M}_{\odot}$,
* Number of haloes: $6784818$.

**Reader:**

We provide a python module: `voids_read.py`, with two functions to read the void catalogues:
* `read_sph` :: for the spherical void finder,
* `read_pop` :: for the popcorn void finder.

We also provide a script: `abundance.py`, to make a quick plot of the radii distribution to check the void identification.
The same scripts are also written in R: `read_voids.r` and `abundance.r`.
