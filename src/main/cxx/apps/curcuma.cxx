/// @file curcuma.cxx
/// Especia for intergalactic metal and non-damped H I, He I, II lines.
/// @author Ralf Quast
/// @date 2023
/// @copyright MIT License
#include <exception>
#include <iostream>

#include "../core/model.h"
#include "../core/runner.h"

using namespace std;


/// Especia to analyse intergalactic metal lines.
///
/// @param argc The number of command line arguments.
/// @param argv[0] The program name.
/// @param argv[1] The random seed.
/// @param argv[2] The parent number.
/// @param argv[3] The population size.
/// @param argv[4] The initial global step size.
/// @param argv[5] The accuracy goal.
/// @param argv[6] The stop generation number.
/// @param argv[7] The trace modulus.
/// @return an exit code
///
/// @remark Usage: curcuma {seed} {parents} {population} {step} {accuracy} {stop} {trace} < {model file} [> {result file}]
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
