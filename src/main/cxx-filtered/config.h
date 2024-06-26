/// @file config.h
/// Configuration constants.
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#ifndef ESPECIA_CONFIG_H
#define ESPECIA_CONFIG_H

#include <string>

namespace especia
{

/// The project name.
const std::string project_name = "@PROJECT_NAME@"; // NOLINT

/// The project title.
const std::string project_title = "@PROJECT_TITLE@"; // NOLINT

/// The project version.
const std::string project_version = "@PROJECT_VERSION@"; // NOLINT

/// The project version control identifier.
const std::string project_revision = "@PROJECT_REVISION@"; // NOLINT

/// The project digital object identifier (DOI).
const std::string project_doi = "@PROJECT_DOI@"; // NOLINT

/// A project digital object identifier (DOI) HTML snippet.
const std::string project_doi_html = "@PROJECT_DOI_HTML@"; // NOLINT

/// The project URL.
const std::string project_url = "@PROJECT_URL@"; // NOLINT

/// The project name, version, and revision.
const std::string project_long_name
    = project_revision.empty() ? "@PROJECT_NAME@-@PROJECT_VERSION@" : "@PROJECT_NAME@-@PROJECT_VERSION@ @PROJECT_REVISION@"; // NOLINT

/// The composite name of the operating system the project is compiled for.
const std::string system_name = "@CMAKE_SYSTEM@"; // NOLINT

/// The vendor of the C++ compiler used to compile the project.
const std::string cxx_compiler = "@CMAKE_CXX_COMPILER_ID@"; // NOLINT

/// The version of the C++ compiler used to compile the project.
const std::string cxx_compiler_version
    = "@CMAKE_CXX_COMPILER_VERSION@"; // NOLINT
} // namespace especia

#endif // ESPECIA_CONFIG_H
