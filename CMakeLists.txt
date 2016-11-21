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

cmake_minimum_required(VERSION 3.6)
project(especia LANGUAGES CXX)
set(PROJECT_VERSION_MAJOR 2017)
set(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_AFFIX snapshot)
set(CMAKE_CONFIGURATION_TYPES Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
    set(VECLIB ${LAPACK_LIBRARIES})
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(MAIN ${CMAKE_SOURCE_DIR}/src/main)
set(TEST ${CMAKE_SOURCE_DIR}/src/test)

configure_file (
        ${MAIN}/cxx/config.h.in
        ${MAIN}/cxx/config.h)

configure_file (
        ${TEST}/html/example.html.in
        ${TEST}/html/example.html)

set(COMMON_SOURCES
        ${MAIN}/cxx/config.h
        ${MAIN}/cxx/model.h
        ${MAIN}/cxx/mtwister.h
        ${MAIN}/cxx/optimize.h
        ${MAIN}/cxx/profiles.cxx
        ${MAIN}/cxx/profiles.h
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

add_executable(xtractcom ${MAIN}/cxx/xtractcom.cxx)
add_executable(xtractdat ${MAIN}/cxx/xtractdat.cxx)
add_executable(xtractlog ${MAIN}/cxx/xtractlog.cxx)
add_executable(xtractmod ${MAIN}/cxx/xtractmod.cxx)
add_executable(xtractmes ${MAIN}/cxx/xtractmes.cxx)

add_custom_target(check)
add_dependencies(check check_example)

add_custom_target(example
        ./xtractmod < ${TEST}/html/example.html | `./xtractcom < ${TEST}/html/example.html` > example.html)
add_dependencies(example especid xtractcom xtractmod)

add_custom_target(check_example
        diff ${TEST}/html/example.html example.html
        COMMENT "Error statements issued by this target and caused by small differences in numerics are normal.")
add_dependencies(check_example example)

add_custom_target(clean_example
        rm -f example.html)

install(TARGETS especia especid especis especiv xtractcom xtractdat xtractlog xtractmod xtractmes
        RUNTIME DESTINATION ~/bin)


add_executable(airtovac EXCLUDE_FROM_ALL ${MAIN}/cxx/util/airtovac.cxx)
add_executable(cumulate EXCLUDE_FROM_ALL ${MAIN}/cxx/util/cumulate.cxx)
add_executable(helicorr EXCLUDE_FROM_ALL ${MAIN}/cxx/util/helicorr.cxx)
add_executable(threecol EXCLUDE_FROM_ALL ${MAIN}/cxx/util/threecol.cxx)
add_executable(vactoair EXCLUDE_FROM_ALL ${MAIN}/cxx/util/vactoair.cxx)