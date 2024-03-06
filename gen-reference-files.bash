#!/bin/bash

SOURCE=`pwd`
cd /tmp
mkdir -p build-hwode
cd build-hwode/
cmake -S $SOURCE
make

RUSSELL="/home/dorival/01-Code/rust/russell/russell_ode/data"

cd examples

./dop853_vdp_debug   > "$RUSSELL/fortran_dop853_van_der_pol_debug.txt"
./dop853_vdp         > "$RUSSELL/fortran_dop853_van_der_pol.txt"
./dopri5_aren_debug  > "$RUSSELL/fortran_dopri5_arenstorf_debug.txt"
./dopri5_aren        > "$RUSSELL/fortran_dopri5_arenstorf.txt"
./dopri5_eq1         > "$RUSSELL/fortran_dopri5_hairer_wanner_eq1.txt"
./dopri5_rob         > "$RUSSELL/fortran_dopri5_robertson.txt"
./dopri5_vdp_debug   > "$RUSSELL/fortran_dopri5_van_der_pol_debug.txt"
./radau5_amp1t       > "$RUSSELL/fortran_radau5_amplifier1t.txt"
./radau5_amp2t_debug > "$RUSSELL/fortran_radau5_amplifier2t_debug.txt"
./radau5_amp2t       > "$RUSSELL/fortran_radau5_amplifier2t.txt"
./radau5_eq1_debug   > "$RUSSELL/fortran_radau5_hairer_wanner_eq1_debug.txt"
./radau5_eq1         > "$RUSSELL/fortran_radau5_hairer_wanner_eq1.txt"
./radau5_rob_debug   > "$RUSSELL/fortran_radau5_robertson_debug.txt"
./radau5_rob_small_h > "$RUSSELL/fortran_radau5_robertson_small_h.txt"
./radau5_rob         > "$RUSSELL/fortran_radau5_robertson.txt"
./radau5_vdp_debug   > "$RUSSELL/fortran_radau5_van_der_pol_debug.txt"
./radau5_vdp         > "$RUSSELL/fortran_radau5_van_der_pol.txt"
