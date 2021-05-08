//! @file eeet.cxx
//! Utility to merge separated spectral flux and uncertainty data
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include <fstream>

#include "../core/base.h"
#include "../core/dataio.h"
#include "../core/exitcodes.h"

using namespace std;

using especia::natural;
using especia::real;

/// Writes the usage message to an output stream.
///
/// @param os The output stream.
/// @param program_name The program name.
void write_usage_message(ostream &os, const string &program_name) {
    os << "usage: " << program_name << " {flux file} {uncertainty file} [lines to skip] [> {target file}]" << endl;
}

/// Utility to merge separate spectral flux and uncertainty data files.
///
/// @param argc The number of command line arguments supplied.
/// @param argv The command line arguments:
/// @parblock
/// @c argv[0] The program name.
///
/// @c argv[1] The path name of the flux data file.
///
/// @c argv[2] The path name of the flux uncertainty data file.
///
/// @c argv[3] The number of lines to skip at the beginning (optional, default = 0).
/// @endparblock
/// @return an exit code.
///
/// @remark Usage: eeet {flux file} {uncertainty file} [lines to skip] [> {target file}]
int main(int argc, char *argv[]) {
    const string program_name(argv[0]);

    if (argc == 1) {
        write_usage_message(cout, program_name);
        return 0;
    }
    try {
        if (argc != 3 and argc != 4) {
            throw invalid_argument("Error: an invalid number of arguments was supplied");
        }

        natural skip = 0;

        if (argc == 4) {
            skip = especia::convert<natural>(string(argv[2]));
        }

        valarray<real> x;
        valarray<real> y;
        valarray<real> z;

        ifstream fxy(argv[1]);
        ifstream fxz(argv[2]);

        especia::get(fxy, x, y, skip);
        especia::get(fxz, x, z, skip);

        fxy.close();
        fxz.close();

        if (fxy and fxz and y.size() == z.size()) {
            especia::put(cout, x, y, z);
        } else {
            throw runtime_error("Error: an input error occurred");
        }
        return 0;
    } catch (logic_error &e) {
        cerr << e.what() << endl;
        return especia::Exit_Codes::logic_error;
    } catch (runtime_error &e) {
        cerr << e.what() << endl;
        return especia::Exit_Codes::runtime_error;
    } catch (exception &e) {
        cerr << e.what() << endl;
        return especia::Exit_Codes::unspecific_exception;
    }
}
