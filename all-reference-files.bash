#!/bin/bash

SOURCE=`pwd`
cd /tmp
mkdir -p build-hwode
cd build-hwode/
cmake -S $SOURCE
make

RUSSELL="/home/dorival/01-Code/rust/russell/russell_ode/data"

cd examples

./dr_dop853       > "$RUSSELL/fortran_dop853_van_der_pol.txt"

./dr_dopri5_debug > "$RUSSELL/fortran_dopri5_arenstorf_debug.txt"
./dr_dopri5       > "$RUSSELL/fortran_dopri5_arenstorf.txt"

./dr_dopri5_eq1   > "$RUSSELL/fortran_dopri5_hairer_wanner_eq1.txt"

./dr0_radau5_debug > "$RUSSELL/fortran_radau5_hairer_wanner_eq1_debug.txt"
./dr0_radau5       > "$RUSSELL/fortran_radau5_hairer_wanner_eq1.txt"

./dr1_radau5_debug > "$RUSSELL/fortran_radau5_van_der_pol_debug.txt"
./dr1_radau5       > "$RUSSELL/fortran_radau5_van_der_pol.txt"
