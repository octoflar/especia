## @author Ralf Quast
## @date 2023
## @copyright MIT License

cmake_minimum_required(VERSION 3.5.0)
include(src/main/cmake/project.cmake)
include(src/main/cmake/unix.cmake)
include(src/main/cmake/test.cmake)
include(src/main/cmake/veclib.cmake)
include(src/main/cmake/openmp.cmake)

project(especia VERSION 2024.1 LANGUAGES C CXX)
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
configure_file(${MAIN}/scripts-filtered/update-copyright update-copyright @ONLY)
configure_file(${MAIN}/scripts-filtered/update-copyright.sed update-copyright.sed @ONLY)

configure_file(${TEST}/resources/curcuma_test.html resources/curcuma_test.html COPYONLY)
configure_file(${TEST}/resources/oregano_test.html resources/oregano_test.html COPYONLY)
configure_file(${TEST}/resources/mascada_test.html resources/mascada_test.html COPYONLY)
configure_file(${TEST}/resources/azafran_test.html resources/azafran_test.html COPYONLY)
configure_file(${TEST}/resources/doublet_test.html resources/doublet_test.html COPYONLY)
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
        ${MAIN}/cxx/core/optimizer.cxx
        ${MAIN}/cxx/core/optimizer.h
        ${MAIN}/cxx/core/optimize.h
        ${MAIN}/cxx/core/profiles.cxx
        ${MAIN}/cxx/core/profiles.h
        ${MAIN}/cxx/core/readline.cxx
        ${MAIN}/cxx/core/readline.h
        ${MAIN}/cxx/core/random.h
        ${MAIN}/cxx/core/runner.cxx
        ${MAIN}/cxx/core/runner.h
        ${MAIN}/cxx/core/section.cxx
        ${MAIN}/cxx/core/section.h)

add_executable(azafran ${MAIN}/cxx/apps/azafran.cxx ${CORE_SOURCES})
target_link_libraries(azafran ${VECLIB} ${CMAKE_THREAD_LIBS_INIT})
add_executable(curcuma ${MAIN}/cxx/apps/curcuma.cxx ${CORE_SOURCES})
target_link_libraries(curcuma ${VECLIB} ${CMAKE_THREAD_LIBS_INIT})
add_executable(oregano ${MAIN}/cxx/apps/oregano.cxx ${CORE_SOURCES})
target_link_libraries(oregano ${VECLIB} ${CMAKE_THREAD_LIBS_INIT})
add_executable(mascada ${MAIN}/cxx/apps/mascada.cxx ${CORE_SOURCES})
target_link_libraries(mascada ${VECLIB} ${CMAKE_THREAD_LIBS_INIT})

add_executable(ecom ${MAIN}/cxx/apps/ecom.cxx)
add_executable(edat ${MAIN}/cxx/apps/edat.cxx)
add_executable(elog ${MAIN}/cxx/apps/elog.cxx)
add_executable(emod ${MAIN}/cxx/apps/emod.cxx)
add_executable(emes ${MAIN}/cxx/apps/emes.cxx)
add_executable(ezip ${MAIN}/cxx/apps/ezip.cxx
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
        ${MAIN}/cxx/core/optimize.h
        ${MAIN}/cxx/core/optimizer.h
        ${MAIN}/cxx/core/optimizer.cxx
        ${TEST}/cxx/core/optimizer_test.cxx
        ${MAIN}/cxx/core/random.h)
target_link_libraries(optimizer_test ${VECLIB} ${CMAKE_THREAD_LIBS_INIT})
add_unit_test(profiles_test
        ${MAIN}/cxx/core/base.h
        ${MAIN}/cxx/core/integrator.h
        ${MAIN}/cxx/core/profiles.h
        ${MAIN}/cxx/core/profiles.cxx
        ${TEST}/cxx/core/profiles_test.cxx)
add_unit_test(random_test
        ${MAIN}/cxx/core/base.h
        ${MAIN}/cxx/core/random.h
        ${TEST}/cxx/core/random_test.cxx)

add_flavour_test(curcuma_test 159.77 171.89)
add_flavour_test(oregano_test 251.04)
add_flavour_test(mascada_test 251.46)
add_flavour_test(azafran_test 103.91 258.67 213.20 223.51 238.06 196.48)

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


install(TARGETS azafran curcuma oregano mascada ecom edat elog emod emes ezip
        RUNTIME
        DESTINATION bin
        CONFIGURATIONS Release)

install(PROGRAMS ${MAIN}/scripts/erun
        DESTINATION bin
        CONFIGURATIONS Release)

add_custom_target(uninstall
        COMMAND ${XARGS} ${RM} < install_manifest.txt
        VERBATIM)
