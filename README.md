# Fortran codes for ordinary differential equations

This repository includes code from Hairer, Norsett and Wanner books.

Here, we present the CMake toolchain for compiling all code.

Please check out [Prof Ernst Hairer' website](http://www.unige.ch/~hairer/software.html)

References:

```
[1] Hairer, Norsett and Wanner (1993): Solving Ordinary Differential Equations.
    Nonstiff Problems. 2nd edition. Springer Series in Comput. Math., vol. 8.
[2] Hairer and Wanner (1996): Solving Ordinary Differential Equations.
    Stiff and Differential-Algebraic Problems. 2nd edition.
    Springer Series in Comput. Math., vol. 14.
```

## Nonstiff Differential Equations

1. DOPRI5 explicit Runge-Kutta method of order 5(4) for problems y'=f(x,y); with dense output of order 4
2. DOP853 explicit Runge-Kutta method of order 8(5,3) for problems y'=f(x,y); with dense output of order 7
3. ODEX Extrapolation method (GBS) for problems y'=f(x,y); with dense output
4. ODEX2 Extrapolation method (Stoermers rule) for second order differential equations y''=f(x,y); with dense output
5. RETARD Extension of the code DOPRI5 to delay differential equations y'(t)=f(t,y(t),y(t-a),...)

## Stiff Differential Equations and Differential-Algebraic Problems

6. RADAU5 implicit Runge-Kutta method of order 5 (Radau IIA) for problems of the form My'=f(x,y) with
   possibly singular matrix M; with dense output (collocation solution).
7. RADAU implicit Runge-Kutta method (Radau IIA) of variable order (switches automatically between
   orders 5, 9, and 13) for problems of the form My'=f(x,y) with possibly singular matrix M; For the
   choices IWORK(11)=3 and IWORK(12)=3, the code is mathematically equivalent to RADAU5 (in general a
   little bit slower than RADAU5).
8. RODAS Rosenbrock method of order 4(3), for problems of the form My'=f(x,y) with possibly singular
   matrix M; with dense output; algebraic order conditions are considered
9. SEULEX Extrapolation method based on linearly implicit Euler for problems of the form My'=f(x,y)
   with possibly singular matrix M; with dense output

### Oldies

10. SDIRK4 a diagonally-implicit Runge-Kutta method of order 4 for problems of the form My'=f(x,y) with possibly singular matrix M; with dense output

## Usage

The code is easily compiled using CMake.

On a Debian-based Linux system, type:

```bash
./all.bash
```

The examples can be run from `/tmp/build-hwode/examples` and `/tmp/build-hwode/tests`.
