# New in 2017.1

The main programs and executable files were renamed and some utility 
programs were added. Ancillary files were revised. Version 2017.1 
implements some new features and enhancements:

* Implementation of the extended pseudo-Voigt approximation (Ida, Ando &
Toraya 2000).
* Replacement of the residual sum of squares with the *cost function*,
as used in inverse problem theory. The cost function simply is half the
residual sum of squares.
* More accurate and faster computation of parameter uncertainties and
the covariance matrix.
* Comment lines in spectrum data files are permitted. Comment lines must
start with a hash mark `#` in the first column.
* The command line typed to invoke Especia is echoed to the standard
output stream. It appears in the result HTML file within a comment block.
The command line can be extracted from the result file with the `ecom`
utility.
* The model definition put in is echoed to standard output. It appears in the
result HTML file within a comment block. The model definition can be extracted
from the result file by means of the `emod` utility.
* The project uses [CMake](https://cmake.org) instead of
[make](https://www.gnu.org/software/make/) to configure and manage the build
and installation processes.

# New in 2016.1

Version 2016.1 is the initial public release. The original source from 2006
was adapted to C++11 coding standards and some compiler warnings eliminated.

