cmake_minimum_required(VERSION 3.6)
project(especia VERSION 2017.1 LANGUAGES CXX)

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

set(MAIN_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/main)
set(TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/test)

configure_file (
        ${TEST_SOURCE_DIR}/resources/example.html.in
        ${TEST_SOURCE_DIR}/resources/example.html
)

set(COMMON_SOURCES
        ${MAIN_SOURCE_DIR}/cxx/profiles.cxx
        ${MAIN_SOURCE_DIR}/cxx/readline.cxx
        ${MAIN_SOURCE_DIR}/cxx/section.cxx
        ${MAIN_SOURCE_DIR}/cxx/symeig.cxx)

add_executable(especia ${MAIN_SOURCE_DIR}/cxx/especia.cxx
        ${COMMON_SOURCES})
target_link_libraries(especia
        ${VECLIB})
add_executable(especid ${MAIN_SOURCE_DIR}/cxx/especid.cxx
        ${COMMON_SOURCES})
target_link_libraries(especid
        ${VECLIB})
add_executable(especis ${MAIN_SOURCE_DIR}/cxx/especis.cxx
        ${COMMON_SOURCES})
target_link_libraries(especis
        ${VECLIB})
add_executable(especiv ${MAIN_SOURCE_DIR}/cxx/especiv.cxx
        ${COMMON_SOURCES})
target_link_libraries(especiv
        ${VECLIB})

add_executable(xtractcom ${MAIN_SOURCE_DIR}/cxx/xtractcom.cxx)
add_executable(xtractdat ${MAIN_SOURCE_DIR}/cxx/xtractdat.cxx)
add_executable(xtractlog ${MAIN_SOURCE_DIR}/cxx/xtractlog.cxx)
add_executable(xtractmod ${MAIN_SOURCE_DIR}/cxx/xtractmod.cxx)
add_executable(xtractmes ${MAIN_SOURCE_DIR}/cxx/xtractmes.cxx)

add_custom_target(example
        ./xtractmod < ${TEST_SOURCE_DIR}/resources/example.html | `./xtractcom < ${TEST_SOURCE_DIR}/resources/example.html` > example.html)
add_custom_target(check
        diff ${TEST_SOURCE_DIR}/resources/example.html example.html)
add_dependencies(check example)

install(TARGETS especia especid especis especiv xtractcom xtractdat xtractlog xtractmod xtractmes
        RUNTIME DESTINATION ~/bin)


add_executable(airtovac EXCLUDE_FROM_ALL ${MAIN_SOURCE_DIR}/cxx/util/airtovac.cxx)
add_executable(cumulate EXCLUDE_FROM_ALL ${MAIN_SOURCE_DIR}/cxx/util/cumulate.cxx)
add_executable(helicorr EXCLUDE_FROM_ALL ${MAIN_SOURCE_DIR}/cxx/util/helicorr.cxx)
add_executable(threecol EXCLUDE_FROM_ALL ${MAIN_SOURCE_DIR}/cxx/util/threecol.cxx)
add_executable(vactoair EXCLUDE_FROM_ALL ${MAIN_SOURCE_DIR}/cxx/util/vactoair.cxx)
