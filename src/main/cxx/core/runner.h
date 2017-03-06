/// @file runner.h
/// The model runner.
/// Copyright (c) 2016 Ralf Quast
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in all
/// copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
/// SOFTWARE.
#ifndef ESPECIA_RUNNER_H
#define ESPECIA_RUNNER_H

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config.h"
#include "optimizer.h"

namespace especia {

    /**
     * Carries out an optimization run.
     */
    class Runner {
    public:

        /**
         * Constructs a new runner from the command line arguments supplied.
         *
         * @param argc The number of command line arguments.
         * @param argv The command line arguments:
         * @parblock
         * @c argv[0] The program name.
         *
         * @c argv[1] The random seed.
         *
         * @c argv[2] The parent number.
         *
         * @c argv[3] The population size.
         *
         * @c argv[4] The initial global step size.
         *
         * @c argv[5] The accuracy goal.
         *
         * @c argv[6] The stop generation number.
         *
         * @c argv[7] The trace modulus.
         * @endparblock
         */
        Runner(int argc, char *argv[]);

        /**
         * Destructor.
         */
        ~Runner();

        /**
         * Returns the command line arguments.
         *
         * @return the command line arguments.
         */
        const std::vector<std::string> &get_args() const {
            return args;
        }

        /**
         * Returns the number of command line arguments.
         *
         * @return the number of command line arguments.
         */
        size_t get_arg_count() const {
            return args.size();
        }

        /**
         * Returns the program name.
         *
         * @return the program name.
         */
        std::string get_program_name() const {
            return std::string(args[0]);
        }

        /**
         * Parses the accuracy goal.
         *
         * @return the accuracy goal.
         */
        Real_t parse_accuracy_goal() const throw(std::invalid_argument) {
            return convert<Real_t>(args[5]);
        }

        /**
         * Parses the initial global step size.
         *
         * @return the inititial global step size.
         */
        Real_t parse_global_step_size() const throw(std::invalid_argument) {
            return convert<Real_t>(args[4]);
        }

        /**
         * Parses the parent number.
         *
         * @return the parent number.
         */
        Nnum_t parse_parent_number() const throw(std::invalid_argument) {
            return convert<Nnum_t>(args[2]);
        }

        /**
         * Parses the population size.
         *
         * @return the population size.
         */
        Nnum_t parse_population_size() const throw(std::invalid_argument) {
            return convert<Nnum_t>(args[3]);
        }

        /**
         * Parses the random seed.
         *
         * @return the random seed.
         */
        Word_t parse_random_seed() const throw(std::invalid_argument) {
            return convert<Word_t>(args[1]);
        }

        /**
         * Parses the stop generation.
         *
         * @return the stop generation.
         */
        Lnum_t parse_stop_generation() const throw(std::invalid_argument) {
            return convert<Lnum_t>(args[6]);
        }

        /**
         * Parses the trace modulus.
         *
         * @return the trace modulus.
         */
        Nnum_t parse_trace_modulus() const throw(std::invalid_argument) {
            return convert<Nnum_t>(args[7]);
        }

        /**
         * Runs the model supplied as argument.
         *
         * @tparam M The model type.
         *
         * @return an exit code
         * @throw invalid_argument when an invalid argument was supplied.
         * @throw runtime_error when a runtime error occurred.
         */
        template<class M>
        int run() throw(std::invalid_argument, std::runtime_error) {
            using std::cin;
            using std::cout;
            using std::endl;
            using std::invalid_argument;
            using std::runtime_error;

            if (get_arg_count() == 1) {
                write_usage_message(cout);
                return 0;
            }
            if (get_arg_count() != 8) {
                throw invalid_argument(
                        "especia::Runner::run() Error: an invalid number of arguments was supplied");
            }

            write_command_line(cout);

            const Word_t random_seed = parse_random_seed();
            const Nnum_t parent_number = parse_parent_number();
            const Nnum_t population_size = parse_population_size();
            const Real_t global_step_size = parse_global_step_size();
            const Real_t accuracy_goal = parse_accuracy_goal();
            const Lnum_t stop_generation = parse_stop_generation();
            const Nnum_t trace_modulus = parse_trace_modulus();

            M model;
            model.get(cin, cout);

            if (cin.fail()) {
                throw runtime_error(
                        "especia::Runner::run() Error: an error occurred while reading the model definition");
            }
            if (not cin.eof()) {
                throw runtime_error(
                        "especia::Runner::run() Error: an error occurred while reading the model definition");
            }

            Optimizer optimizer = Optimizer::Builder().
                    with_problem_dimension(model.get_parameter_count()).
                    with_parent_number(parent_number).
                    with_population_size(population_size).
                    with_accuracy_goal(accuracy_goal).
                    with_stop_generation(stop_generation).
                    with_random_seed(random_seed).
                    build();

            cout << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">" << endl;
            cout << "<html>" << endl;

            cout << "<!--" << endl;
            cout << "<log>" << endl;

            Optimizer::Result result = optimizer.minimize(model,
                                                          model.get_initial_parameter_values(),
                                                          model.get_initial_local_step_sizes(),
                                                          global_step_size,
                                                          model.get_constraint(),
                                                          Tracer<>(cout, trace_modulus));

            cout << "</log>" << endl;
            cout << "-->" << endl;

            write_result_messages(cout, result);

            cout << "</html>" << endl;

            model.set(&result.get_parameter_values()[0], &result.get_parameter_uncertainties()[0]);
            model.put(cout);

            if (result.is_optimized()) {
                return 0;
            } else {
                return 1;
            }
        }

    private:
        /**
         * Traces optimizer state information to an output stream.
         *
         * @tparam T The number type.
         */
        template<class T = Real_t>
        class Tracer {
        public:
            /**
             * Constructor.
             *
             * @param[in] output_stream The output stream.
             * @param[in] modulus The trace modulus.
             * @param[in] precision The precision of numeric output.
             * @param[in] width The width of the numeric output fields.
             */
            Tracer(std::ostream &output_stream, Nnum_t modulus, Nnum_t precision = 4, Nnum_t width = 12)
                    : os(output_stream), m(modulus), p(precision), w(width) {
            }

            /**
             * Destructor.
             */
            ~Tracer() {
            }

            /**
             * Tests if tracing is enabled.
             *
             * @param[in] g The generation number.
             * @return @true if tracing is enabled, otherwise @c false.
             */
            Bool_t is_enabled(Lnum_t g) const {
                return m > 0 and g % m == 0;
            }

            /**
             * Traces state information to an output stream..
             *
             * @param[in] g The generation number.
             * @param[in] y The value of the objective function.
             * @param[in] min_step The minimum step size.
             * @param[in] max_step The maximum step size.
             */
            void trace(Lnum_t g, T y, T min_step, T max_step) const {
                using std::endl;
                using std::ios_base;
                using std::setw;

                const ios_base::fmtflags fmt = os.flags();

                os.setf(ios_base::fmtflags());
                os.setf(ios_base::scientific, ios_base::floatfield);
                os.setf(ios_base::right, ios_base::adjustfield);
                os.precision(p);

                os << setw(8) << g;
                os << setw(w) << y;
                os << setw(w) << min_step;
                os << setw(w) << max_step;
                os << endl;

                os.flags(fmt);
            }

        private:
            std::ostream &os;
            const Nnum_t m;
            const Nnum_t p;
            const Nnum_t w;
        };

        void write_command_line(std::ostream &os) const;

        void write_result_messages(std::ostream &os, const Optimizer::Result &result) const;

        void write_usage_message(std::ostream &os) const;

        /**
         * The command line arguments.
         */
        std::vector<std::string> args;
    };

}

#endif //ESPECIA_RUNNER_H
