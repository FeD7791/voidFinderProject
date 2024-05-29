# VoidFinder Toolkit
The Void Finder Toolkit (VFT is a Python software package that integrates various publicly available algo-
rithms to characterize voids and represent them in a standardized and comparable manner, irrespective
of their geometry. The implementation includes a common interface inspired by the popular Scikit-Learn
architecture [Buitinck et al., 2013] that enables users, using the same set of tracers, to determine whether a
tracer belongs to a void as defined by a selected algorithm, facilitating a more direct comparison of results.

## Description
The code integrates three cosmic void detection algorithms, implemented in C/C++ to harness low-level
execution speed, with a common interface written in Python that offers users expressiveness and flexibil-
ity for interactive analysis.

### Finders

The current status of VFT integrates 3 public algorithms:
• Zobov [Neyrinck, 2008]: ZOBOV works in analogy with a watershed method with water filling basins
in a density field. It looks for voids as density minima with surrounding depressions and requires no
free parameters. Each Void grows in density starting from a local minimum up to a link density where
particles start falling into a deeper minimum.
• Spherical [Ruiz et al., 2015]: This method searches regions of low density in a Voronoi tessellation.
For each minimum density region the algorithm then grows a sphere around each candidate until the
average density inside reaches a specific threshold.
• Popcorn [Paz et al., 2023]: The algorithm targets low-density regions by adding layers on spherical void
shapes. Each layer strategically places seeds that expand while maintaining density. Only the best seed
merges, and a refined process ensures full coverage. This continues until small spheres can’t be added,
capturing the entire void effectively.

**INPUT**

In order to simplify the representation of astronomical input data, we utilize the concept of a Box object
type to encapsulate all the complexity related to the tracers. So, in our context, a box is an object that
contains a set of tracers along with their properties, serving as input for the different void finder methods.
The input dataset must align with the structure defined by this object. Each box consists of three coordi-
nate components, three tracer velocity components, and their corresponding masses. Currently, this input
structure is designed for data sourced from a catalog of tracers (galaxies or dark matter) obtained from
cosmological simulations.

**METHOD**
The project internally involves a collection of modules, objects, classes, and functions that collectively ful-
fill the objective of providing an interactive and computationally efficient environment for the analysis of
Voids using different algorithms.
To simplify the reader’s understanding, we have decided to present the computational architecture we
propose by explaining how the data “moves” from the tracers to obtaining the final result.
1. A finder is created by configuring the selected low-level algorithm with the desired parameters (each
algorithm has a set of configurable parameters, which may include).
2. The user provides the Box containing all the tracer information to the algorithm.
3. The tracer preprocesses the box and transforms it into a suitable structure to feed the low-level algo-
rithm.
4. The actual void search algorithm is executed using the parameters provided in step 1, and it is fed with
the preprocessed tracer data from step 3.
5. The data returned by the algorithm is post-processed and remapped back to the tracers in a matrix
where each row represents a tracer and each column represents a void. If a tracer belongs to a void, the
corresponding cell will have a value of 1, and 0 otherwise.
6. A VoidedBox object is created, containing the original box, the matrix from step 5, and any interesting
data resulting from the algorithm in step 4.
7. The VoidedBox is returned to the user.






