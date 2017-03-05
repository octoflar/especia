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

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"

namespace especia {

    /**
     * Runs a model.
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
         * @c argv[7] The trace interval.
         * @endparblock
         */
        Runner(int argc, char *argv[]);

        /**
         * Destructor.
         */
        ~Runner();

        /**
         * Returns the accuracy goal.
         *
         * @return the accuracy goal.
         */
        double get_accuracy_goal() const {
            return static_cast<double>(std::atof(args[5].c_str()));
        }

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
        unsigned long get_arg_count() const {
            return args.size();
        }

        /**
         * Returns the initial global step size.
         *
         * @return the inititial global step size.
         */
        double get_global_step_size() const {
            return static_cast<double>(std::atof(args[4].c_str()));
        }

        /**
         * Returns the parent number.
         *
         * @return the parent number.
         */
        unsigned int get_parent_number() const {
            return static_cast<unsigned>(std::atoi(args[2].c_str()));
        }

        /**
         * Returns the population size.
         *
         * @return the population size.
         */
        unsigned int get_population_size() const {
            return static_cast<unsigned>(std::atoi(args[3].c_str()));
        }

        /**
         * Returns the program name.
         *
         * @return the program name.
         */
        std::string get_program_name() const {
            return args[0];
        }

        /**
         * Returns the random seed.
         * @return the random seed.
         */
        unsigned long get_random_seed() const {
            return static_cast<unsigned long>(std::atol(args[1].c_str()));
        }

        /**
         * Returns the stop generation.
         *
         * @return the stop generation.
         */
        unsigned long get_stop_generation() const {
            return static_cast<unsigned long>(std::atol(args[6].c_str()));
        }

        /**
         * Returns the trace interval.
         *
         * @return the trace interval.
         */
        unsigned int get_trace_interval() const {
            return static_cast<unsigned>(std::atoi(args[7].c_str()));
        }

        /**
         * Runs the model supplied as argument.
         *
         * @tparam M The model type.
         *
         * @param model The model.
         * @return an exit code: 0 = OK, 1 = model is not optimized, 2-4 = input error
         *
         * A usage message is witten to standard output, if no command line arguments (excluding
         * the program name) were supplied. In this case zero is returned.
         */
        template<class M>
        int run(M &model) {
            if (get_arg_count() == 1) {
                write_usage_message(std::cout);

                return 0;
            }
            if (get_arg_count() != 8) {
                write_usage_message(std::cout);

                return 2;
            }

            write_command_line(std::cout);

            const unsigned long seed = get_random_seed();
            const unsigned int parent_number = get_parent_number();
            const unsigned int population_size = get_population_size();
            const double step_size = get_global_step_size();
            const double accuracy_goal = get_accuracy_goal();
            const unsigned long stop = get_stop_generation();
            const unsigned int trace = get_trace_interval();

            model.get(std::cin, std::cout);

            if (std::cin.fail()) {
                return 3;
            }
            if (not std::cin.eof()) {
                return 4;
            }

            const bool optimized = model.optimize(parent_number,
                                                  population_size,
                                                  step_size,
                                                  accuracy_goal,
                                                  stop,
                                                  trace,
                                                  seed, std::cout);
            model.put(std::cout);

            if (not optimized) {
                return 1;
            }

            return 0;
        }

    private:
        void write_command_line(std::ostream &os) const {
            os << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">" << std::endl;
            os << "<html>" << std::endl;
            os << "<!--" << std::endl;
            os << "<command>" << std::endl;

            for (int i = 0; i < args.size(); ++i) {
                os << " " << args[i];
            }

            os << std::endl;
            os << "</command>" << std::endl;
            os << "-->" << std::endl;
            os << "</html>" << std::endl;
        }

        void write_usage_message(std::ostream &os) const {
            os << PROJECT_LONG_NAME << " " << DOI << std::endl;
            os << "usage: " << get_program_name() << ": "
               << "SEED PARENTS POPULATION INISTEP ACCURACY STOPGEN TRACE < ISTREAM > OSTREAM"
               << std::endl;
        }

        /**
         * The command line arguments.
         */
        std::vector<std::string> args;
    };

}

#endif //ESPECIA_RUNNER_H
