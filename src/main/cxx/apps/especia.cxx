/// @file especia.cxx
/// Especia for inferring the variation of the fine-structure constant.
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include <exception>
#include <iostream>

#include "../core/model.h"
#include "../core/runner.h"

using namespace std;


/// Flavor of Especia to infer the variation of the fine-structure constant.
///
/// @param argc The number of command line arguments.
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
/// @remark Usage: especia {seed} {parents} {population} {step} {accuracy} {stop} {trace} < {model file} [> {result file}]
///
/// @attention A usage message is written to standard output, if no command line arguments (excluding the
/// program name) are supplied. In this case the returned exit code is zero.
int main(int argc, char *argv[]) {
    typedef especia::Model<especia::Many_Multiplet> Model;

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
