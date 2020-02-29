# Copyright (c) 2020 Ralf Quast
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.2.2)
include(src/main/cmake/project.cmake)
include(src/main/cmake/unix.cmake)
include(src/main/cmake/test.cmake)
include(src/main/cmake/veclib_required.cmake)
include(src/main/cmake/openmp_optional.cmake)

project(especia VERSION 2020.1 LANGUAGES C CXX)
project_version_tag(snapshot)
project_title("Evolutionary spectrum inversion and analysis (Especia)")
project_doi("https://doi.org/10.5281/zenodo.785071")
project_doi_html("<a href=\\\"https://doi.org/10.5281/zenodo.785071\\\">doi.org/10.5281/zenodo.785071</a>")
project_url("https://github.com/octoflar/especia")
project_install_prefix($ENV{HOME})

set(CMAKE_CONFIGURATION_TYPES Release Debug)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES config.h)

veclib_required()
openmp_optional()

set(MAIN ${CMAKE_SOURCE_DIR}/src/main)
set(TEST ${CMAKE_SOURCE_DIR}/src/test)
include_directories(${CMAKE_BINARY_DIR})

configure_file(${MAIN}/cxx-filtered/config.h config.h @ONLY)
configure_file(${MAIN}/scripts/erun erun COPYONLY)

configure_file(${TEST}/resources/especid_test.html resources/especid_test.html COPYONLY)
configure_file(${TEST}/resources/especiv_test.html resources/especiv_test.html COPYONLY)
configure_file(${TEST}/resources/especix_test.html resources/especix_test.html COPYONLY)
configure_file(${TEST}/resources/especia_test.html resources/especia_test.html COPYONLY)
configure_file(${TEST}/resources/supersampling_test.html resources/supersampling_test.html COPYONLY)
configure_file(${TEST}/resources/spec000.dat resources/spec000.dat COPYONLY)
configure_file(${TEST}/resources/stis121.dat resources/stis121.dat COPYONLY)
configure_file(${TEST}/resources/uves346.dat resources/uves346.dat COPYONLY)
configure_file(${TEST}/resources/uves580.dat resources/uves580.dat COPYONLY)
configure_file(${TEST}/resources/virt156.dat resources/virt156.dat COPYONLY)
configure_file(${TEST}/resources/virt166.dat resources/virt166.dat COPYONLY)


set(CORE_SOURCES
        ${MAIN}/cxx/core/base.h
        config.h
        ${MAIN}/cxx/core/decompose.cxx
        ${MAIN}/cxx/core/decompose.h
        ${MAIN}/cxx/core/deviates.h
        ${MAIN}/cxx/core/exitcodes.h
        ${MAIN}/cxx/core/integrator.h
        ${MAIN}/cxx/core/model.h
        ${MAIN}/cxx/core/mtwister.h
        ${MAIN}/cxx/core/optimizer.cxx
        ${MAIN}/cxx/core/optimizer.h
        ${MAIN}/cxx/core/optimize.h
        ${MAIN}/cxx/core/profiles.cxx
        ${MAIN}/cxx/core/profiles.h
        ${MAIN}/cxx/core/readline.cxx
        ${MAIN}/cxx/core/readline.h
        ${MAIN}/cxx/core/runner.cxx
        ${MAIN}/cxx/core/runner.h
        ${MAIN}/cxx/core/section.cxx
        ${MAIN}/cxx/core/section.h)

add_executable(especia ${MAIN}/cxx/apps/especia.cxx ${CORE_SOURCES})
target_link_libraries(especia ${VECLIB})
add_executable(especid ${MAIN}/cxx/apps/especid.cxx ${CORE_SOURCES})
target_link_libraries(especid ${VECLIB})
add_executable(especiv ${MAIN}/cxx/apps/especiv.cxx ${CORE_SOURCES})
target_link_libraries(especiv ${VECLIB})
add_executable(especix ${MAIN}/cxx/apps/especix.cxx ${CORE_SOURCES})
target_link_libraries(especix ${VECLIB})
add_executable(ecom ${MAIN}/cxx/apps/ecom.cxx)
add_executable(edat ${MAIN}/cxx/apps/edat.cxx)
add_executable(elog ${MAIN}/cxx/apps/elog.cxx)
add_executable(emod ${MAIN}/cxx/apps/emod.cxx)
add_executable(emes ${MAIN}/cxx/apps/emes.cxx)
add_executable(eeet ${MAIN}/cxx/apps/eeet.cxx
        ${MAIN}/cxx/core/dataio.cxx
        ${MAIN}/cxx/core/dataio.h
        ${MAIN}/cxx/core/exitcodes.h)

add_unit_test(decompose_test
        ${MAIN}/cxx/core/base.h
        ${MAIN}/cxx/core/decompose.h
        ${MAIN}/cxx/core/decompose.cxx
        ${TEST}/cxx/core/decompose_test.cxx)
target_link_libraries(decompose_test ${VECLIB})
add_unit_test(integrator_test
        ${MAIN}/cxx/core/base.h
        ${MAIN}/cxx/core/integrator.h
        ${TEST}/cxx/core/integrator_test.cxx)
add_unit_test(optimizer_test
        ${MAIN}/cxx/core/base.h
        ${MAIN}/cxx/core/decompose.h
        ${MAIN}/cxx/core/decompose.cxx
        ${MAIN}/cxx/core/deviates.h
        ${MAIN}/cxx/core/mtwister.h
        ${MAIN}/cxx/core/optimize.h
        ${MAIN}/cxx/core/optimizer.h
        ${MAIN}/cxx/core/optimizer.cxx
        ${TEST}/cxx/core/optimizer_test.cxx)
target_link_libraries(optimizer_test ${VECLIB})
add_unit_test(profiles_test
        ${MAIN}/cxx/core/base.h
        ${MAIN}/cxx/core/integrator.h
        ${MAIN}/cxx/core/profiles.h
        ${MAIN}/cxx/core/profiles.cxx
        ${TEST}/cxx/core/profiles_test.cxx)

add_sanity_test(supersampling_test 13.89)
add_performance_test(especid_test 159.77 171.89)
add_performance_test(especiv_test 251.04)
add_performance_test(especix_test 251.46)
add_performance_test(especia_test 103.91 258.67 213.20 223.51 238.06 196.48)

add_executable(airtovac EXCLUDE_FROM_ALL ${MAIN}/cxx/util/airtovac.cxx
        ${MAIN}/cxx/core/base.h
        ${MAIN}/cxx/core/dataio.cxx
        ${MAIN}/cxx/core/dataio.h
        ${MAIN}/cxx/core/equations.cxx
        ${MAIN}/cxx/core/equations.h
        ${MAIN}/cxx/core/exitcodes.h)
add_executable(helicorr EXCLUDE_FROM_ALL ${MAIN}/cxx/util/helicorr.cxx
        ${MAIN}/cxx/core/base.h
        ${MAIN}/cxx/core/dataio.cxx
        ${MAIN}/cxx/core/dataio.h
        ${MAIN}/cxx/core/exitcodes.h)
add_executable(vactoair EXCLUDE_FROM_ALL ${MAIN}/cxx/util/vactoair.cxx
        ${MAIN}/cxx/core/base.h
        ${MAIN}/cxx/core/dataio.cxx
        ${MAIN}/cxx/core/dataio.h
        ${MAIN}/cxx/core/equations.cxx
        ${MAIN}/cxx/core/equations.h
        ${MAIN}/cxx/core/exitcodes.h)

add_custom_target(util)
add_dependencies(util airtovac helicorr vactoair)

add_custom_target(tar
        COMMAND ${TAR} -f ${PROJECT_NAME}-${PROJECT_VERSION}-${PROJECT_VERSION_TAG}.tar.gz -C ${CMAKE_SOURCE_DIR} --exclude .git --exclude .gitattributes --exclude .gitignore --exclude .travis.yml --exclude _config.yml -X ${CMAKE_SOURCE_DIR}/.gitignore .
        BYPRODUCTS ${PROJECT_NAME}-${PROJECT_VERSION}-${PROJECT_VERSION_TAG}.tar.gz
        VERBATIM)


install(TARGETS especia especid especiv especix ecom edat elog emod emes eeet
        RUNTIME
        DESTINATION bin
        CONFIGURATIONS Release)

install(PROGRAMS ${MAIN}/scripts/erun
        DESTINATION bin
        CONFIGURATIONS Release)

add_custom_target(uninstall
        COMMAND ${XARGS} ${RM} < install_manifest.txt
        VERBATIM)
