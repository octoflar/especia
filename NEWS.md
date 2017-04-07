# New in 2017.1

The main programs were renamed and a few utility programs were added.
The source code, documentation, and ancillary files were revised.
Version 2017.1 provides new features and enhancements:

* The equivalent width of each line is calculated and listed in the
result table.
* Equations and equation solvers to convert photon wavelength in vacuum
to photon wavelength in standard air (and vice versa) were added to the
core API (Edlén 1953, 1966; Birch & Downs 1994).
* Especia uses standard C++-11 multi-threading when the compiler does
not support Open Multiprocessing directives.
* The spectral resolution of the instrument is expressed in units of
10<sup>3</sup> which improves the scaling of the optimization problem.
* The variation of the fine-structure constant Δα/α is calculated in
units of 10<sup>-6</sup> which improves the scaling of the problem.
* The extended pseudo-Voigt approximation (Ida, Ando &
Toraya 2000) was implemented.
* The residual sum of squares was replaced with the *cost function*,
termed in inverse problem theory. The cost function simply is half the
residual sum of squares.
* The computation of model parameter uncertainties and covariance was
improved in terms of accuracy and speed.
* Comment lines in spectrum data files are permitted now. Comment lines
must start with a hash mark `#`, an exclamation mark `!`, or a percent
sign `%` in the first column.
* The command line typed to invoke Especia is included in the result HTML
file within a comment block. This command line can be extracted from the
result file by means of the new `ecom` utility.
* The model definition is included in the result HTML file within a comment
block. The model definition can be extracted from the result file by means
of the new `emod` utility.
* GNU [make](https://www.gnu.org/software/make/) was replaced with the
cross-platform [CMake](https://cmake.org) to make the build and installation
processes simpler.

# New in 2016.1

Version 2016.1 is the initial public release. The original source from 2006
was adapted to C++11 coding standards and some compiler warnings eliminated.

