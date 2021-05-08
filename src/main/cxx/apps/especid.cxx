//! @file especid.cxx
//! Especia for intergalactic metal and non-damped H I, He I, II lines.
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include <exception>
#include <iostream>

#include "../core/model.h"
#include "../core/runner.h"

using namespace std;


/// Flavor of Especia to analyse intergalactic metal lines.
///
/// @param argc The number of command line arguments.
/// @param argv The command line arguments:
/// @parblock
/// @c argv[0] The program name.
///
/// @c argv[1] The random seed.
///
/// @c argv[2] The parent number.
///
/// @c argv[3] The population size.
///
/// @c argv[4] The initial global step size.
///
/// @c argv[5] The accuracy goal.
///
/// @c argv[6] The stop generation number.
///
/// @c argv[7] The trace modulus.
/// @endparblock
/// @return an exit code
///
/// @remark Usage: especid {seed} {parents} {population} {step} {accuracy} {stop} {trace} < {model file} [> {result file}]
///
/// @attention A usage message is written to standard output, if no command line arguments (excluding the
/// program name) are supplied. In this case the returned exit code is zero.
int main(int argc, char *argv[]) {
    typedef especia::Model<especia::Intergalactic_Doppler> Model;

    try {
        return especia::Runner(argc, argv).run<Model>();
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
