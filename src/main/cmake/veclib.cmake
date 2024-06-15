## @author Ralf Quast
## @date 2024
## @copyright MIT License

macro(veclib_required)
    find_package(LAPACK REQUIRED)
    if (LAPACK_FOUND)
        set(VECLIB ${LAPACK_LIBRARIES} CACHE STRING
                "The libraries for doing linear algebra")
    endif ()
endmacro()
