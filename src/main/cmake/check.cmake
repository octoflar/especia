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

function(add_check NAME EXECUTABLE RESOURCE EXPECTED)
    get_filename_component(BASENAME ${RESOURCE} NAME)
    add_custom_command(OUTPUT ${BASENAME}
            COMMAND ./emod < ${RESOURCE} | `./ecom < ${RESOURCE}` > ${BASENAME})
    add_custom_target(${NAME}
            DEPENDS ${BASENAME}
            COMMENT "The check is passed, if the next lines issue the number 1.")
    foreach (VALUE ${EXPECTED} ${ARGN})
        add_custom_command(TARGET ${NAME} PRE_BUILD
                COMMAND ${GREP} --count --fixed-strings '<td><strong>${VALUE}</strong></td>' ${BASENAME})
    endforeach ()
    add_dependencies(${NAME} ecom emod ${EXECUTABLE})
    add_dependencies(check ${NAME})
endfunction()

add_custom_target(check)
