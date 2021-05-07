//! @file config.h
//! Configuration constants.

// Copyright (c) 2021. Ralf Quast
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#ifndef ESPECIA_CONFIG_H
#define ESPECIA_CONFIG_H

#include <string>

namespace especia {

    /**
     * The project name.
     */
    const std::string project_name = "@PROJECT_NAME@"; // NOLINT

    /**
     * The project title.
     */
    const std::string project_title = "@PROJECT_TITLE@"; // NOLINT

    /**
     * The project version.
     */
    const std::string project_version = "@PROJECT_VERSION@"; // NOLINT

    /**
     * The project version control tag.
     */
    const std::string project_version_tag = "@PROJECT_VERSION_TAG@"; // NOLINT

    /**
     * The project digital object identifier (DOI).
     */
    const std::string project_doi = "@PROJECT_DOI@"; // NOLINT

    /**
     * A project digital object identifier (DOI) HTML snippet.
     */
    const std::string project_doi_html = "@PROJECT_DOI_HTML@"; // NOLINT

    /**
     * The project URL.
     */
    const std::string project_url = "@PROJECT_URL@"; // NOLINT

    /**
     * The project name and version identifier.
     */
    const std::string project_long_name = "@PROJECT_NAME@-@PROJECT_VERSION@ @PROJECT_VERSION_TAG@"; // NOLINT

    /**
     * The composite name of the operating system the project is compiled for.
     */
    const std::string system_name = "@CMAKE_SYSTEM@"; // NOLINT

    /**
     * The vendor of the C++ compiler used to compile the project.
     */
    const std::string cxx_compiler = "@CMAKE_CXX_COMPILER_ID@"; // NOLINT

    /**
     * The version of the C++ compiler used to compile the project.
     */
    const std::string cxx_compiler_version = "@CMAKE_CXX_COMPILER_VERSION@"; // NOLINT
}

#endif // ESPECIA_CONFIG_H
