## @author Ralf Quast
## @date 2024
## @copyright MIT License

macro(openmp_optional)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} $<$<CONFIG:Release>:${OpenMP_C_FLAGS}>")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $<$<CONFIG:Release>:${OpenMP_CXX_FLAGS}>")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} $<$<CONFIG:Release>:${OpenMP_EXE_LINKER_FLAGS}>")
    else ()
        find_package(Threads REQUIRED)
    endif ()
endmacro()
