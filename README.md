# Synopsis

The Evolutionary spectrum inversion and analysis (Especia) file set
provides ISO C++ code for the inverse modelling and analysis of intergalactic and
interstellar absorption line regions in QSO spectra.

Due to the use of evolution strategies with covariance matrix adaption, the inverse
modelling procedure is highly competitive and capable of calculating the optimal
spectral decomposition without requiring any particular initialisation.

Further highlights are the modelling and optimisation of the background continuum,
and a very accurate semi-analytic convolution of the absorption term with the
instrumental profile.

The algorithms are explained in
[Quast et al. (2005)](http://dx.doi.org/10.1051/0004-6361:20041601).
A provided example case uses an artificial spectrum synthesised on basis of
data described and analysed by
[Quast et al. (2002)](http://dx.doi.org/10.1051/0004-6361:20020342).


# Getting started

This software enables you to analyse spectroscopic data. Though it has been developed
for the analysis of astrophysical spectra, it is applicable to spectroscopy in general.
Read the two articles listed below to find out whether this software is of interest to
you. If you already have some understanding of these matters, consult the
[especia wiki](https://github.com/octoflar/especia/wiki) for operating instructions.

Building this software requires [CMake](https://cmake.org) and a compiler that implements
the ISO/IEC 14882:2011 norm, also known as C++11. To build and test this software
`cd` into the project root directory and type:

    mkdir cmake-build
    cd cmake-build
    cmake ..
    make check

Typing `make install` will build all software and move the executable files into your
`$HOME/bin` directory. In case of problems consult the
[build instructions](https://github.com/octoflar/especia/wiki/Build-instructions).


# Versioning

Release versions YYYY.N are numbered by the year of the release follwowed by a
single-digit number, which enumerates the release within the release year. For
example, version 2016.1 denotes the first release of the year 2016.


# Further reading

Quast, Ralf; Baade, Robert; Reimers, Dieter (2005): *Evolution strategies applied to the problem of line profile decomposition in QSO spectra.*
Astronomy and Astrophysics 431 (3) 1167.
doi: [10.1051/0004-6361:20041601](http://dx.doi.org/10.1051/0004-6361:20041601).

Quast, Ralf; Baade, Robert; Reimers, Dieter (2002): *Fine-structure diagnostics of neutral carbon toward HE 0515-4414.*
Astronomy and Astrophysics 386 (3) 796.
doi: [10.1051/0004-6361:20020342](http://dx.doi.org/10.1051/0004-6361:20020342).
