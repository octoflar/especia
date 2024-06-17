## @author Ralf Quast
## @date 2024
## @copyright MIT License

function(add_feature_test NAME)
    add_test(NAME ${NAME} COMMAND ./erun resources/${NAME}.html ${NAME}.html ${ARGN} CONFIGURATIONS Release)
    set_tests_properties(${NAME} PROPERTIES LABELS feature TIMEOUT 10)
    add_custom_target(${NAME} ctest --output-on-failure --tests-regex ${NAME})
endfunction()

add_custom_target(featuretests ctest --output-on-failure --label-regex feature)

function(add_flavour_test NAME)
    add_test(NAME ${NAME} COMMAND ./erun resources/${NAME}.html ${NAME}.html ${ARGN} CONFIGURATIONS Release)
    set_tests_properties(${NAME} PROPERTIES LABELS flavour TIMEOUT 900)
    add_custom_target(${NAME} ctest --verbose --tests-regex ${NAME})
endfunction()

add_custom_target(flavourtests ctest --verbose --label-regex flavour)

function(add_unit_test NAME)
    add_executable(${NAME}
            ${ARGN}
            ${TEST}/cxx/unittest.h)
    target_compile_options(${NAME} PRIVATE $<$<CONFIG:Debug>:--coverage>)
    target_link_options(${NAME} PRIVATE $<$<CONFIG:Debug>:--coverage>)
    add_test(NAME ${NAME} COMMAND ${NAME})
    set_tests_properties(${NAME} PROPERTIES LABELS unit TIMEOUT 10)
endfunction()

add_custom_target(unittests ctest --output-on-failure --label-regex unit)

enable_testing()
