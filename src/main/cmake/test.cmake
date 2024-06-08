## @author Ralf Quast
## @date 2021
## @copyright MIT License

function(add_flavour_test NAME EXPECTED_VALUES)
    add_test(NAME ${NAME} COMMAND ./erun resources/${NAME}.html ${NAME}.html ${EXPECTED_VALUES} ${ARGN})
    set_tests_properties(${NAME} PROPERTIES LABELS flavour TIMEOUT 900)
    add_custom_target(${NAME} ctest --verbose --tests-regex ${NAME})
endfunction()

add_custom_target(flavourtests ctest --verbose --label-regex flavour)

function(add_unit_test NAME)
    add_executable(${NAME}
            ${ARGN}
            ${TEST}/cxx/unittest.h)
    target_compile_options(${NAME} PRIVATE --coverage)
    target_link_options(${NAME} PRIVATE --coverage)
    add_test(NAME ${NAME} COMMAND ${NAME})
    set_tests_properties(${NAME} PROPERTIES LABELS unit TIMEOUT 10)
endfunction()

add_custom_target(unittests ctest --output-on-failure --label-regex unit)

enable_testing()
