cmake_minimum_required(VERSION 3.5)
project(hwode)

option(OPTIMIZED "optimized (-O2)" OFF)
option(VERBOSE "verbose" OFF)

if(OPTIMIZED)
	add_definitions(-O2)
else(OPTIMIZED)
    add_definitions(-g)
endif(OPTIMIZED)

if(VERBOSE)
	set(CMAKE_VERBOSE_MAKEFILE TRUE)
endif(VERBOSE)

enable_language(Fortran)

include(FindLAPACK)

set(LIBS    ${LIBS}   ${LAPACK_LIBRARIES})
set(LFLAGS "${LFLAGS} ${LAPACK_LINKER_FLAGS}")

file(GLOB all *.f)

add_library(hwode STATIC
    hinit.f
    dc_decsol.f
    dc_lapack.f
    decsol.f
    dop853.f
    dopri5.f
    odex2.f
    odex.f
    radau5.f
    radau.f
    retard.f
    rodas.f
    sdirk4.f
    seulex.f
)
target_link_libraries(hwode ${LIBS})
set_target_properties(hwode PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")

set(LFLAGS "-L.")
set(LIBS ${LIBS} hwode)

subdirs(examples tests)
