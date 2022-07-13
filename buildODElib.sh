#!/bin/bash

# Create object files in ode/lib/
cd ./ode/lib/

gcc -I ../include/ -c ../src/euler_solver.c
gcc -I ../include/ -c ../src/rk2_solver.c
gcc -I ../include/ -c ../src/rk4_solver.c


# Create library
ar -r libODE.a euler_solver.o rk2_solver.o rk4_solver.o