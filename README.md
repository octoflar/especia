## Synopsis

The Evolutionary spectrum inversion and analysis (Especia) package provides
C++ code for the inverse modelling and analysis of intergalactic and interstellar
absorption line regions in QSO spectra.

Includes the modelling and optimisation of the background continuum by a linear
combination of Legendre polynomials. Implements an accurate semi-analytic convolution
of the absorption term with the instrumental profile. Also includes IO for model
definition (plain text) and result (HTML) files.

The code is further explained in
[Quast et al. (2005)](http://dx.doi.org/10.1051/0004-6361:20041601).
The provided example and test case uses an artificial spectrum synthesised on basis
of data described and analysed by
[Quast et al. (2002)](http://dx.doi.org/10.1051/0004-6361:20020342).

## Getting Started

To build test and install the software type

    make
    make test
    make install

Then study the example files and read the literature cited below to find out whether
this software may be useful for you.

## Versioning

Release versions YYYY.N are numbered by the year of the release follwowed by an
single-digit number, which enumerates the release within the release year. For
example, version 2016.1 denotes the first release within the year 2016.

## Authors

The author of the Especia code is 

**Ralf Quast**

affiliated with (1998 - 2006)

*Universität Hamburg, Hamburger Sternwarte, 21029 Hamburg, Germany*

## License

This project is licensed under [The MIT License (MIT)](http://opensource.org/licenses/MIT).

## Acknowledgements

The use of this software in the works listed below (and a few others not listed) is gratefully
acknowledged.

1. Fine-structure diagnostics of neutral carbon toward  HE 0515-4414.
   R. Quast, R. Baade and D. Reimers.
   Astronomy and Astrophysics 386 (3) 796 (2002).
   DOI: [10.1051/0004-6361:20020342](http://dx.doi.org/10.1051/0004-6361:20020342)
2. Probing the variability of the fine-structure constant with the VLT/UVES.
   R. Quast, D. Reimers and S. A. Levshakov.
   Astronomy and Astrophysics 415 (2) L7 (2004).
   DOI: [10.1051/0004-6361:20040013](http://dx.doi.org/10.1051/0004-6361:20040013)
3. Evolution strategies applied to the problem of line profile decomposition in QSO spectra.
   R. Quast, R. Baade and D. Reimers.
   Astronomy and Astrophysics 431 (3) 1167 (2005).
   DOI: [10.1051/0004-6361:20041601](http://dx.doi.org/10.1051/0004-6361:20041601)
4. HE 0515–4414: an unusual sub-damped Ly α system revisite.d
   R. Quast, D. Reimers and R. Baade.
   Astronomy and Astrophysics 477 (2) 443 (2008).
   DOI: [10.1051/0004-6361:20054773](http://dx.doi.org/10.1051/0004-6361:20054773)
5. Robust limit on a varying proton-to-electron mass ratio from a single H2 system.
   M. Wendt and P. Molaro.
   Astronomy & Astrophysics 526 A96 (2011).
   DOI: [10.1051/0004-6361/201014835](http://dx.doi.org/10.1051/0004-6361/201014835)
6. The precision of line position measurements of unresolved quasar absorption lines and its influence on the search for variations of fundamental constants.
   N. Prause and D. Reimers.
   Astronomy & Astrophysics 555 A88 (2013).
   DOI: [10.1051/0004-6361/201118373](http://dx.doi.org/10.1051/0004-6361/201118373)
