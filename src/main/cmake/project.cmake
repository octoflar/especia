# Copyright (c) 2017 Ralf Quast
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

function(project_url URL)
    set(PROJECT_URL ${URL} PARENT_SCOPE)
endfunction()

function(project_tag TAG)
    if (${TAG} STREQUAL snapshot)
        find_program(GIT git)
        if (GIT)
            execute_process(COMMAND ${GIT} rev-parse --short HEAD
                    OUTPUT_VARIABLE ID
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE STATUS
                    ERROR_QUIET)
            if (${STATUS} EQUAL 0)
                set(PROJECT_TAG ${ID} PARENT_SCOPE)
            endif ()
        endif ()
    else ()
        set(PROJECT_TAG ${TAG} PARENT_SCOPE)
    endif ()
endfunction()
