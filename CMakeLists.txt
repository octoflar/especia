# Copyright (c) 2016 Ralf Quast
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
include(src/main/cmake/project_version_tag.cmake)
include(src/main/cmake/veclib_required.cmake)
include(src/main/cmake/openmp_optional.cmake)

project(especia VERSION 2017.1 LANGUAGES C CXX)
project_version_tag(snapshot)
set(PROJECT_DOI https://dx.doi.org/10.6084/m9.figshare.4167999.v14)

set(CMAKE_CONFIGURATION_TYPES Release)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

veclib_required()
openmp_optional()

set(MAIN ${CMAKE_SOURCE_DIR}/src/main)
set(TEST ${CMAKE_SOURCE_DIR}/src/test)

configure_file(${MAIN}/cxx/config.h.in ${MAIN}/cxx/config.h)
configure_file(${MAIN}/resources/synth.html.in resources/synth.html)
configure_file(${MAIN}/resources/synth.md resources/synth.md COPYONLY)
configure_file(${MAIN}/resources/synth1561.dat resources/synth1561.dat COPYONLY)
configure_file(${MAIN}/resources/synth1657.dat resources/synth1657.dat COPYONLY)

set(COMMON_SOURCES
        ${MAIN}/cxx/base.h
        ${MAIN}/cxx/config.h
        ${MAIN}/cxx/model.h
        ${MAIN}/cxx/mtwister.h
        ${MAIN}/cxx/optimize.h
        ${MAIN}/cxx/profiles.cxx
        ${MAIN}/cxx/profiles.h
        ${MAIN}/cxx/randev.h
        ${MAIN}/cxx/readline.cxx
        ${MAIN}/cxx/readline.h
        ${MAIN}/cxx/section.cxx
        ${MAIN}/cxx/section.h
        ${MAIN}/cxx/symeig.cxx
        ${MAIN}/cxx/symeig.h)

add_executable(especia ${MAIN}/cxx/especia.cxx ${COMMON_SOURCES})
target_link_libraries(especia ${VECLIB})
add_executable(especid ${MAIN}/cxx/especid.cxx ${COMMON_SOURCES})
target_link_libraries(especid ${VECLIB})
add_executable(especis ${MAIN}/cxx/especis.cxx ${COMMON_SOURCES})
target_link_libraries(especis ${VECLIB})
add_executable(especiv ${MAIN}/cxx/especiv.cxx ${COMMON_SOURCES})
target_link_libraries(especiv ${VECLIB})
add_executable(ecom ${MAIN}/cxx/ecom.cxx)
add_executable(edat ${MAIN}/cxx/edat.cxx)
add_executable(elog ${MAIN}/cxx/elog.cxx)
add_executable(emod ${MAIN}/cxx/emod.cxx)
add_executable(emsg ${MAIN}/cxx/emsg.cxx)

add_custom_target(check
        COMMAND ./emod < resources/synth.html | `./ecom < resources/synth.html` > check.html
        COMMAND grep --count --fixed-strings '<td><strong>159.77</strong></td>' check.html
        COMMAND grep --count --fixed-strings '<td><strong>171.89</strong></td>' check.html
        BYPRODUCTS check.html
        COMMENT "The check is passed, if the next two lines issue the number 1.")
add_dependencies(check ecom emod especid)

install(TARGETS especia especid especis especiv ecom edat elog emod emsg
        RUNTIME
        DESTINATION $ENV{HOME}/bin
        CONFIGURATIONS Release)

add_executable(airtovac EXCLUDE_FROM_ALL ${MAIN}/cxx/util/airtovac.cxx ${MAIN}/cxx/base.h)
add_executable(barycorr EXCLUDE_FROM_ALL ${MAIN}/cxx/util/barycorr.cxx ${MAIN}/cxx/base.h)
add_executable(cumulate EXCLUDE_FROM_ALL ${MAIN}/cxx/util/cumulate.cxx)
add_executable(threecol EXCLUDE_FROM_ALL ${MAIN}/cxx/util/threecol.cxx)
add_executable(vactoair EXCLUDE_FROM_ALL ${MAIN}/cxx/util/vactoair.cxx)

add_custom_target(util)
add_dependencies(util airtovac barycorr cumulate threecol vactoair)
