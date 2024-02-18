set(EXAMPLES
    dr0_radau5
    dr1_radau5
    dr2_radau5
    dr_dop853
    dr_dopri5
    dr_odex2
    dr_odex
    dr_radau
    dr_retard
    dr_rodas
    dr_seulex
)

foreach(var ${EXAMPLES})
    add_executable       (${var} "${var}.f")
    target_link_libraries(${var} ${LIBS})
    set_target_properties(${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
endforeach(var)
