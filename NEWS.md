# New in 2017.1

Version 2017.1 is a maintenance release. The main programs and executable
files were renamed and some utility programs were added. Ancillary files
were revised. Build properties for the GCC and Intel compilers were added.
A few new features and enhancements were implemented:

* More accurate and faster computation of parameter uncertainties and error
covariance
* Comment lines in spectrum data files are permitted. Comment lines must
start with a hash mark `#` in the first column.
* The command line typed to invoke `especi*` is echoed to the standard
output stream. It appears in the result HTML file within a comment block.
The command line can be extracted from the result file with the `xtractcom`
utility.
* The model definition put in is echoed to standard output. It appears in the
result HTML file within a comment block. The model definition can be extracted
from the result file by means of the `xtractmod` utility.
* The project uses [CMake](https://cmake.org) instead of
[make](https://www.gnu.org/software/make/) to configure and manage the build
processes.

# New in 2016.1

Version 2016.1 is the initial public release. The original sources from 2006
were adapted to C++11 coding standards and some compiler warnings eleminated.

