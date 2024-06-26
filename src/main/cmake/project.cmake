## @author Ralf Quast
## @date 2024
## @copyright MIT License

function(project_doi DOI)
    set(PROJECT_DOI ${DOI} PARENT_SCOPE)
endfunction()

function(project_doi_html HTML)
    set(PROJECT_DOI_HTML ${HTML} PARENT_SCOPE)
endfunction()

function(project_title TITLE)
    set(PROJECT_TITLE ${TITLE} PARENT_SCOPE)
endfunction()

function(project_url URL)
    set(PROJECT_URL ${URL} PARENT_SCOPE)
endfunction()

function(project_rev REV)
    if (${REV} STREQUAL GIT)
        find_program(GIT git)
        if (GIT)
            execute_process(COMMAND ${GIT} rev-parse --short HEAD
                    OUTPUT_VARIABLE REV
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE STATUS
                    ERROR_QUIET)
            if (${STATUS} EQUAL 0)
                set(PROJECT_REVISION ${REV} PARENT_SCOPE)
            else ()
                set(PROJECT_REVISION "" PARENT_SCOPE)
            endif ()
        endif ()
    else ()
        set(PROJECT_REVISION ${REV} PARENT_SCOPE)
    endif ()
endfunction()

function(project_install_prefix PREFIX)
    if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "${PREFIX}" CACHE PATH
                "This directory is prepended onto all install directories"
                FORCE)
    endif ()
endfunction()
