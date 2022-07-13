#! /bin/bash

# Rebuild exponential decay executable
gcc -I ode/include/ decay.c -o decay.exe -L ode/lib/ -lODE

# Generate data, gets put in `decay_data` directory
./decay.exe 0.01 0.02 0.04 0.05 0.1 0.2 0.25 0.5 1.0