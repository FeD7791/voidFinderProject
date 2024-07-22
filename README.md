# VoidFinder Toolkit
The Void Finder Toolkit VFTK is a Python software package that integrates various publicly available algorithms to characterize voids and represent them in a standardized and comparable manner, irrespective of their geometry. The implementation includes a common interface inspired by the popular Scikit-Learn architecture [Buitinck et al., 2013].

## Description
The code integrates three cosmic void detection algorithms, implemented in C/C++ to harness low-level execution speed, with a common interface written in Python that offers users expressiveness and flexibility for interactive analysis.

### Finders

The current status of VFT integrates 3 public algorithms:
• Zobov [Neyrinck, 2008]: ZOBOV works in analogy with a watershed method with water filling basins in a density field. It looks for voids as density minima with surrounding depressions and requires no free parameters. Each Void grows in density starting from a local minimum up to a link density where particles start falling into a deeper minimum.
• Spherical [Ruiz et al., 2015]: This method searches regions of low density in a Voronoi tessellation. For each minimum density region the algorithm then grows a sphere around each candidate until the average density inside reaches a specific threshold.
• Popcorn [Paz et al., 2023]: The algorithm targets low-density regions by adding layers on spherical void shapes. Each layer strategically places seeds that expand while maintaining density. Only the best seed merges, and a refined process ensures full coverage. This continues until small spheres can’t be added, capturing the entire void effectively.

## dev Installation

Clone this repo and then inside the local directory execute

```bash
$ git clone https://github.com/FeD7791/voidFinderProject.git
$ cd voidFinderProject
$ pip install -e .
```
### Install Zobov
GCC is needed
```bash
$ cd voidFinderProject
$ cd .voidfindertk/zobov/src
$ make
```







