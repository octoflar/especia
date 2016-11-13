## General information

The **E**volutionary **spec**trum *i*nversion and *a*nalysis (Especia) file set
provides C++ code for the inverse modelling and analysis of intergalactic and
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

## Getting Started

To build, test, and install this software type

    make
    make test
    make install

Then study the example files and read the articles listed below to find out whether
this software could be of interest to you. You may also want to consult the
[especia wiki](https://github.com/octoflar/especia/wiki).

## Versioning

Release versions YYYY.N are numbered by the year of the release follwowed by a
single-digit number, which enumerates the release within the release year. For
example, version 2016.1 denotes the first release of the year 2016.

## Further reading

Quast, Ralf; Baade, Robert; Reimers, Dieter (2005): *Evolution strategies applied to the problem of line profile decomposition in QSO spectra.*
Astronomy and Astrophysics 431 (3) 1167.
[DOI: 10.1051/0004-6361:20041601](http://dx.doi.org/10.1051/0004-6361:20041601).

Quast, Ralf; Baade, Robert; Reimers, Dieter (2002): *Fine-structure diagnostics of neutral carbon toward HE 0515-4414.*
Astronomy and Astrophysics 386 (3) 796.
[DOI: 10.1051/0004-6361:20020342](http://dx.doi.org/10.1051/0004-6361:20020342).
