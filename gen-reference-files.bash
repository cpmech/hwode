#!/bin/bash

SOURCE=`pwd`
cd /tmp
mkdir -p build-hwode
cd build-hwode/
cmake -S $SOURCE
make

RUSSELL="/home/dorival/01-Code/rust/russell/russell_ode/data"

cd examples

./dop853_vdp        > "$RUSSELL/fortran_dop853_van_der_pol.txt"
./dopri5_aren_debug > "$RUSSELL/fortran_dopri5_arenstorf_debug.txt"
./dopri5_aren       > "$RUSSELL/fortran_dopri5_arenstorf.txt"
./dopri5_eq1        > "$RUSSELL/fortran_dopri5_hairer_wanner_eq1.txt"
./radau5_amp_debug  > "$RUSSELL/fortran_radau5_amplifier_debug.txt"
./radau5_amp        > "$RUSSELL/fortran_radau5_amplifier.txt"
./radau5_eq1_debug  > "$RUSSELL/fortran_radau5_hairer_wanner_eq1_debug.txt"
./radau5_eq1        > "$RUSSELL/fortran_radau5_hairer_wanner_eq1.txt"
./radau5_rob        > "$RUSSELL/fortran_radau5_robertson.txt"
./radau5_vdp_debug  > "$RUSSELL/fortran_radau5_van_der_pol_debug.txt"
./radau5_vdp        > "$RUSSELL/fortran_radau5_van_der_pol.txt"
