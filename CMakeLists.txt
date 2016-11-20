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

set(COMMON_SOURCES
        src/main/cxx/profiles.cxx
        src/main/cxx/readline.cxx
        src/main/cxx/section.cxx
        src/main/cxx/symeig.cxx)

add_executable(especia src/main/cxx/especia.cxx
        ${COMMON_SOURCES})
target_link_libraries(especia
        ${VECLIB})
add_executable(especid src/main/cxx/especid.cxx
        ${COMMON_SOURCES})
target_link_libraries(especid
        ${VECLIB})
add_executable(especis src/main/cxx/especis.cxx
        ${COMMON_SOURCES})
target_link_libraries(especis
        ${VECLIB})
add_executable(especiv src/main/cxx/especiv.cxx
        ${COMMON_SOURCES})
target_link_libraries(especiv
        ${VECLIB})

add_executable(xtractcom src/main/cxx/xtractcom.cxx)
add_executable(xtractdat src/main/cxx/xtractdat.cxx)
add_executable(xtractlog src/main/cxx/xtractlog.cxx)
add_executable(xtractmod src/main/cxx/xtractmod.cxx)
add_executable(xtractmes src/main/cxx/xtractmes.cxx)

add_custom_target(example
        ./xtractmod < src/test/resources/example.html | `./xtractcom < src/test/resources/example.html` > example.html)
add_custom_target(check
        diff src/test/resources/example.html example.html)
add_dependencies(check example)

install(TARGETS especia especid especis especiv xtractcom xtractdat xtractlog xtractmod xtractmes
        RUNTIME DESTINATION ~/bin)


add_executable(airtovac EXCLUDE_FROM_ALL src/main/cxx/util/airtovac.cxx)
add_executable(cumulate EXCLUDE_FROM_ALL src/main/cxx/util/cumulate.cxx)
add_executable(helicorr EXCLUDE_FROM_ALL src/main/cxx/util/helicorr.cxx)
add_executable(threecol EXCLUDE_FROM_ALL src/main/cxx/util/threecol.cxx)
add_executable(vactoair EXCLUDE_FROM_ALL src/main/cxx/util/vactoair.cxx)
