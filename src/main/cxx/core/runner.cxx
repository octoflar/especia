/// @file runner.h
/// The model runner.
/// Copyright (c) 2017 Ralf Quast
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
#include "runner.h"

especia::Runner::Runner(int argc, char *argv[]) {
    using std::string;

    for (int i = 0; i < argc; ++i) {
        args.emplace_back(argv[i]);
    }
}

especia::Runner::~Runner() = default;

void especia::Runner::write_command_line(std::ostream &os) const {
    using std::endl;

    os << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">" << endl;
    os << "<html>" << endl;
    os << "<!--" << endl;
    os << "<command>" << endl;

    for (const auto &arg : args) {
        os << " " << arg;
    }

    os << std::endl;
    os << "</command>" << endl;
    os << "-->" << endl;
    os << "</html>" << endl;
}

void especia::Runner::write_result_messages(std::ostream &os, const Optimizer::Result &result) const {
    using std::cout;
    using std::endl;

    os << "<!--" << endl;
    os << "<message>" << endl;

    if (result.is_optimized()) {
        os << "especia::Runner::run() Message: optimization completed successfully"
           << endl;
    } else {
        os << "especia::Runner::run() Warning: optimization stopped at generation "
           << result.get_generation_number()
           << endl;
    }
    if (result.is_underflow()) {
        os << "especia::Runner::run() Warning: optimization stopped because of an "
           << "underflow of the mutation variance"
           << endl;
    }

    os << "</message>" << endl;
    os << "-->" << endl;
}

void especia::Runner::write_usage_message(std::ostream &os) const {
    using std::endl;

    os << project_long_name << " " << project_doi << endl;
    os << "usage: " << get_program_name() << ": "
       << "{seed} {parents} {population} {step} {accuracy} {stop} {trace} < {model file} [> {result file}]"
       << endl;
}
