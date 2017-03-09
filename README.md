# Synopsis

The Evolutionary spectrum inversion and analysis (Especia) file set
provides ISO C++ code for the inverse modelling and analysis of intergalactic and
interstellar absorption line regions in QSO spectra.

Due to the use of evolution strategies with covariance matrix adaption, the inverse
modelling procedure is highly competitive and capable of calculating the optimal
spectral decomposition without requiring any particular initialisation.

Additional highlights are the joint modelling and optimisation of the background
continuum, and an accurate semi-analytic convolution of the absorption term with
the instrumental line spread function.

The method is explained in detail by
[Quast et al. (2005)](http://dx.doi.org/10.1051/0004-6361:20041601).


# Getting started [![Build Status](https://travis-ci.org/octoflar/especia.svg?branch=master)](https://travis-ci.org/octoflar/especia)

Especia enables you to analyse spectroscopic data. Though it has been developed
for the analysis of astrophysical spectra, it is, in principle, applicable to
spectroscopy in general. Read the notes and articles listed below to find out
whether this software is of interest to you. If you already have some understanding
of these matters, you may want to consult the 
[especia wiki](https://github.com/octoflar/especia/wiki).

Building Especia requires [CMake](https://cmake.org) and a compiler that implements
[C++11](https://en.wikipedia.org/wiki/C%2B%2B11). To build and test `cd` into the
project root directory and type:

    mkdir cmake-build-release
    cd cmake-build-release
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make all test

Typing `make install` will complete the build and move the executable files into your
`$HOME/bin` directory. For further information or in case of problems consult the
[build instructions](https://github.com/octoflar/especia/wiki/Build-instructions).


# Versions

Release versions YYYY.N are numbered by the year of the release follwowed by a
single-digit number, which enumerates the release within the release year. For
example, version 2016.1 denotes the first release of the year 2016.

The current version is 2016.1.


# Further reading

Quast, Ralf (2014): *Covariance matrix adaption in evolution strategies.* figshare.
doi: [10.6084/m9.figshare.994249.v1](https://doi.org/10.6084/m9.figshare.994249.v1).

Quast, Ralf (2014): *Evolution strategies applied to the problem of line profile decomposition in QSO spectra.* figshare.
doi: [10.6084/m9.figshare.994250.v1](https://doi.org/10.6084/m9.figshare.994250.v1).

Quast, Ralf; Baade, Robert; Reimers, Dieter (2005): *Evolution strategies applied to the problem of line profile decomposition in QSO spectra.*
Astronomy and Astrophysics 431 (3) 1167.
doi: [10.1051/0004-6361:20041601](http://dx.doi.org/10.1051/0004-6361:20041601).
