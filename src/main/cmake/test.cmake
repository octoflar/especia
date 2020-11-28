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

function(add_performance_test NAME EXPECTED_VALUES)
    add_test(NAME ${NAME} COMMAND ./erun resources/${NAME}.html ${NAME}.html ${EXPECTED_VALUES} ${ARGN})
    set_tests_properties(${NAME} PROPERTIES LABELS performance TIMEOUT 3600)
    add_custom_target(${NAME} ctest --verbose --tests-regex ${NAME})
endfunction()

add_custom_target(performancetests ctest --verbose --label-regex performance)

function(add_integration_test NAME EXPECTED_VALUES)
    add_test(NAME ${NAME} COMMAND ./erun resources/${NAME}.html ${NAME}.html ${EXPECTED_VALUES} ${ARGN})
    set_tests_properties(${NAME} PROPERTIES LABELS sanity TIMEOUT 300)
    add_custom_target(${NAME} ctest --verbose --tests-regex ${NAME})
endfunction()

add_custom_target(integrationtests ctest --verbose --label-regex sanity)

function(add_unit_test NAME)
    add_executable(${NAME}
            ${ARGN}
            ${TEST}/cxx/unittest.h)
    add_test(NAME ${NAME} COMMAND ${NAME})
    set_tests_properties(${NAME} PROPERTIES LABELS unit TIMEOUT 10)
endfunction()

add_custom_target(unittests ctest --output-on-failure --label-regex unit)

enable_testing()
