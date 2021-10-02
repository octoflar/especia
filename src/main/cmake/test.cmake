## @author Ralf Quast
## @date 2021
## @copyright MIT License

function(add_property_test NAME EXPECTED_VALUES)
    add_test(NAME ${NAME} COMMAND ./erun resources/${NAME}.html ${NAME}.html ${EXPECTED_VALUES} ${ARGN})
    set_tests_properties(${NAME} PROPERTIES LABELS property TIMEOUT 600)
    add_custom_target(${NAME} ctest --verbose --tests-regex ${NAME})
endfunction()

add_custom_target(propertytests ctest --verbose --label-regex property)

function(add_unit_test NAME)
    add_executable(${NAME}
            ${ARGN}
            ${TEST}/cxx/unittest.h)
    add_test(NAME ${NAME} COMMAND ${NAME})
    set_tests_properties(${NAME} PROPERTIES LABELS unit TIMEOUT 10)
endfunction()

add_custom_target(unittests ctest --output-on-failure --label-regex unit)

enable_testing()
