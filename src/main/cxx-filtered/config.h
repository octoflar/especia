/// @file config.h
/// Configuration constants.
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
#ifndef ESPECIA_CONFIG_H
#define ESPECIA_CONFIG_H

#cmakedefine ESPECIA_PURE_CMAES_VERSION @ESPECIA_PURE_CMAES_VERSION@

#include <string>

namespace especia {

    /**
     * The version of the pure CMA-ES reference implementation. Use '2006' to select the
     * original version.
     */
    const unsigned especia_pure_cmaes_version = @ESPECIA_PURE_CMAES_VERSION@;

    /**
     * The project name.
     */
    const std::string project_name = "@PROJECT_NAME@";

    /**
     * The project version.
     */
    const std::string project_version = "@PROJECT_VERSION@";

    /**
     * The project version control tag.
     */
    const std::string project_version_tag = "@PROJECT_VERSION_TAG@";

    /**
     * The project digital object identifier (DOI).
     */
    const std::string project_doi = "@PROJECT_DOI@";

    /**
     * The project URL.
     */
    const std::string project_url = "@PROJECT_URL@";

    /**
     * The project reference identifier. Either the project DOI, if defined, or the project URL.
     */
    const std::string project_ref = project_doi.empty() ? project_url : project_doi;

    /**
     * The project name and version identifier.
     */
    const std::string project_long_name =
            "@PROJECT_NAME@-@PROJECT_VERSION@ @PROJECT_VERSION_TAG@ (@ESPECIA_PURE_CMAES_VERSION@)";

    /**
     * The composite name of the operating system the project is compiled for.
     */
    const std::string system_name = "@CMAKE_SYSTEM@";

    /**
     * The vendor of the C++ compiler used to compile the project.
     */
    const std::string cxx_compiler = "@CMAKE_CXX_COMPILER_ID@";

    /**
     * The version of the Fortran compiler used to compile the project.
     */
    const std::string cxx_compiler_version = "@CMAKE_CXX_COMPILER_VERSION@";
}

#endif // ESPECIA_CONFIG_H
