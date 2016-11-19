# New in 2017.1

Version 2017.1 is a maintenance release. The main programs and executable
files were renamed and some utility programs were added. Ancillary files
were revised. Build properties for the GCC and Intel compilers were added.
A few new features and enhancements were implemented:

* Comment lines in spectrum data files are permitted. Comment lines must
start with a hash mark `#` in the first column.
* The command line typed to invoke `especi*` runs is echoed to the standard
output stream. It appears in the result HTML file within a comment block.
The command line can be extracted from the result file with the `xtractcom`
utility.
* The model definition read is echoed to standard output. It appears in the
result HTML file within a comment block. The model definition can be extracted
from the result file by means of the `xtractmod` utility.


# New in 2016.1

Version 2016.1 is the initial public release. The source code from 2006 was
adapted to C++11 standards and some compiler warnings were eleminated.
