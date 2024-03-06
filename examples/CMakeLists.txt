set(EXAMPLES
    dop853_vdp_debug
    dop853_vdp
    dopri5_aren_debug
    dopri5_aren
    dopri5_eq1
    dopri5_rob
    dopri5_vdp_debug
    dr_odex
    dr_odex2
    dr_radau
    dr_retard
    dr_rodas
    dr_seulex
    radau5_amp1t
    radau5_amp2t_debug
    radau5_amp2t
    radau5_eq1_debug
    radau5_eq1
    radau5_rob
    radau5_rob_debug
    radau5_rob_small_h
    radau5_vdp_debug
    radau5_vdp
)

foreach(var ${EXAMPLES})
    add_executable       (${var} "${var}.f")
    target_link_libraries(${var} ${LIBS})
    set_target_properties(${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
endforeach(var)
