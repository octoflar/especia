# New in 2024.1

The main programs were renamed and a few utility programs were added.
The source code, documentation, and ancillary files were revised.
In addition, Version 2024.1 provides several new features and enhancements:

* Code test coverage report is generated.
* Especia flavours have been renamed `azafran` for many-multiplet decomposition
  to test the hypothetical variation of the fine-structure constant, `curcuma`
  for Doppler profile decomposition of metal absorption lines, `oregano` for
  Voigt profile decomposition using the pseudo-Voigt approximation, and `pimiento`
  for Voigt profile decomposition using the extended pseudo-Voigt approximation.
* Especia now uses 64-bit instead of 32-bit seeds to initialize the random
  number generator.
* MELG and PCG random number generators are included with the source code, 
  in addition to the Mersenne Twister.
* Especia uses all principal axes of the mutation ellipsoid to estimate parameter 
  uncertainties.
* Especia executes a super-sampled computation of the instrumental convolution 
  to establish accuracy when the spectrum data are not oversampled per se.
* New default values of strategy parameters like recombination weights, 
  and cumulation and adaption rates have been adopted from Nikolaus 
  Hansen's [pure CMA-ES](http://cma.gforge.inria.fr/purecmaes.m) reference 
  implementation.
* The rest equivalent width of each line is calculated and listed in the 
  result table.
* Equations and equation solvers to convert photon wavelength in vacuum 
  to photon wavelength in standard air (and vice versa) are included with 
  the core API (Edlén 1953, 1966; Birch & Downs 1994).
* Standard C++-11 multi-threading is used when the compiler does 
  not support Open Multiprocessing directives.
* The spectral resolution of the instrument is expressed in units of 
  10<sup>3</sup> which improves the scaling of the optimization problem.
* The variation of the fine-structure constant Δα/α is calculated in 
  units of 10<sup>-6</sup> which improves the scaling of the optimization 
  problem.
* The extended pseudo-Voigt approximation (Ida, Ando & Toraya 2000) is 
  implemented.
* The residual sum of squares is replaced with the *cost function*, 
  as used in inverse problem theory. The cost function simply is half the 
  residual sum of squares.
* The computation of model parameter uncertainties and covariance is 
  improved in terms of accuracy and speed.
* Comment lines in spectrum data files are permitted now. Comment lines 
  must start with a hash mark `#`, an exclamation mark `!`, or a percent 
  sign `%` in the first column.
* The new `erun` script runs the command included with a result HTML file.
* The command line typed to invoke Especia is included with the result HTML 
  file within a comment block. The command line can be extracted from the 
  result file by means of the new `ecom` application.
* The model definition is included in the result HTML file within a comment 
  block. The model definition can be extracted from the result file by means 
  of the new `emod` application.
* The cross-platform [CMake](https://cmake.org) replaces GNU make to simplify 
  the build and installation processes.

# New in 2016.1

Version 2016.1 is an initial public pre-release. The source code is almost
identical to the source code of February 2006, except for a few lines of code
that were adapted to C++11 standards to eliminate compiler errors and warnings.
