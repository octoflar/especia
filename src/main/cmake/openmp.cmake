## @author Ralf Quast
## @date 2024
## @copyright MIT License

macro(openmp_optional)
    find_package(OpenMP)
    if ($<$<CONFIG:Release>:OPENMP_FOUND>)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")    
    else ()
        find_package(Threads REQUIRED)
    endif ()
endmacro()
