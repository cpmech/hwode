cmake_minimum_required(VERSION 3.5)
project(hwode)

enable_language(Fortran)

include(FindLAPACK)

set(libs    ${libs}   ${LAPACK_LIBRARIES})
set(lflags "${lflags} ${LAPACK_LINKER_FLAGS}")

file(GLOB all *.f)

add_library(hwode SHARED
    dop853.f
    dopri5.f
    odex2.f
    odex.f
    radau5.f
    radau.f
    retard.f
    rodas.f
    seulex.f
)

#add_executable(calllapack ${all})
#target_link_libraries(calllapack ${libs})
#set_target_properties(calllapack PROPERTIES COMPILE_FLAGS "${flags}" LINK_FLAGS "${lflags}")
