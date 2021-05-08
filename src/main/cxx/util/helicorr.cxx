/// @file helicorr.cxx
/// Utility to apply a heliocentric (or barycentric) velocity correction
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include <exception>
#include <stdexcept>

#include "../core/base.h"
#include "../core/dataio.h"
#include "../core/exitcodes.h"

using namespace std;

using especia::natural;
using especia::real;


/// Writes the usage message to an output stream.
///
/// @param os The ouput stream.
/// @param pname The program name.
void write_usage_message(ostream &os, const string &pname) {
    os << "usage: " << pname << " {velocity (m s-1)} [skip] < {source data file} [> {target data file}]" << endl;
}

/// Utility to apply the heliocentric (or barycentric) velocity correction to
/// spectroscopic data.
///
/// @param argc The number of command line arguments supplied.
/// @param argv The command line arguments:
/// @parblock
/// @c argv[0] The program name.
///
/// @c argv[1] The velocity of the observer relative to the heliocenter (or
/// barycenter) of the solar system (m s-1) projected along the line of sight
/// toward the observed object.
///
/// @c argv[2] The number of lines to skip at the beginning (optional, default = 0).
/// @endparblock
/// @return an exit code.
///
/// @remark Usage: helicorr {velocity (m s-1)} [lines to skip] < {source file} [> {target file}]
int main(int argc, char *argv[]) {
    const string program_name(argv[0]);

    if (argc == 1) {
        write_usage_message(cout, program_name);
        return 0;
    }
    try {
        if (argc != 2 and argc != 3) {
            throw invalid_argument("Error: an invalid number of arguments was supplied");
        }

        const auto v = especia::convert<real>(string(argv[1]));

        natural skip = 0;

        if (argc == 3) {
            skip = especia::convert<natural>(string(argv[2]));
        }

        valarray<real> x;
        valarray<real> y;
        valarray<real> z;

        if (especia::get(cin, x, y, z, skip)) {
            if (v != 0.0) {
                x *= 1.0 + especia::redshift(v);
            }
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
