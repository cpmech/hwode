set(TESTS
    vdp
    rob
    ore
    hir
    e5
    plate
    beam
    ks
    cusp
    brus
    beck
    br2d
)

foreach(var ${TESTS})
    set(bin "${var}-radau5")
    add_executable       (${bin} "${var}_equation.f" "${var}_radau5.f")
    target_link_libraries(${bin} ${LIBS})
    set_target_properties(${bin} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
endforeach(var)

foreach(var ${TESTS})
    set(bin "${var}-seulex")
    add_executable       (${bin} "${var}_equation.f" "${var}_seulex.f")
    target_link_libraries(${bin} ${LIBS})
    set_target_properties(${bin} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
endforeach(var)

