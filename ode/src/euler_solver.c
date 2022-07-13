/* Ravi Chepuri */

#include <stdio.h>
#include "euler_solver.h"


/* SINGLE SPECIES EULER'S METHOD */

void euler_solver(double *x, double *t, double (*dxdt)(double, double, double *),
                    double *dxdt_parameters, double dt, int steps) {
    /* Uses Euler's method to approximate the solution of an arbitrary ODE
                        dx/dt = f(x, t, parameters)
        x: a pointer to an array of size steps to be filled with x values solved 
        t: a pointer to an array of size steps to be filled with t values
        dxdt: a pointer to the function f(x, t, parameters)
        dxdt_parameters: an array that holds parameters for the function f(x, t, parameters)
        dt: time step
        steps: number of steps to iterate Euler's method */

    int i;
    for (i = 0; i < steps - 1; i++) {
        double f = (*dxdt)(*(x+i), *(t+i), dxdt_parameters); // value of f at this timestep
        printf("f: %f\n", f);
        *(x+i+1) = euler_step(*(x+i), f, dt);
        *(t+i+1) = *(t+i) + dt;
    }
}

double euler_step(double x, double f, double dt) {
    return x + f*dt;
}


/* TWO SPECIES EULER'S METHOD */

void euler_solver_twospecies(TwoSpecies *x, double *t,
                    double (*dxAdt)(TwoSpecies, double, double *),
                    double (*dxBdt)(TwoSpecies, double, double *),
                    double *dxdt_parameters, double dt, int steps) {
    /* Uses Euler's method to approximate the solution of an arbitrary two species ODE
                        dx_A/dt = f_A(x_A, x_B, t, parameters)
                        dx_B/dt = f_B(x_A, x_B, t, parameters)
        x: a pointer to a TwoSpecies array of size steps to be filled with x values solved 
        t: a pointer to an array of size steps to be filled with t values
        dxAdt: a pointer to the function f_A(x_A, x_B, t, parameters)
        dxBdt: a pointer to the function f_B(x_A, x_B, t, parameters)
        dxdt_parameters: an array that holds parameters for the function f(x, t, parameters)
        dt: time step
        steps: number of steps to iterate Euler's method */

    int i;
    for (i = 0; i < steps - 1; i++) {
        double f_A = (*dxAdt)(*(x+i), *(t+i), dxdt_parameters); // value of f at this timestep
        double f_B = (*dxBdt)(*(x+i), *(t+i), dxdt_parameters);
        // printf("f_A: %f, f_B: %f\n", f_A, f_B);
        *(x+i+1) = euler_step_twospecies(*(x+i), f_A, f_B, dt);
        *(t+i+1) = *(t+i) + dt;
    }
}

TwoSpecies euler_step_twospecies(TwoSpecies x, double f_A, double f_B, double dt) {
    double new_Na = x.Na + f_A * dt;
    double new_Nb = x.Nb + f_B * dt;
    TwoSpecies new_x = {new_Na, new_Nb};
    return new_x;
}


/* FOUR SPECIES EULER'S METHOD */

void euler_solver_fourspecies(FourSpecies *x, double *t,
                    double (*dx1dt)(FourSpecies, double, double *),
                    double (*dx2dt)(FourSpecies, double, double *),
                    double (*dx3dt)(FourSpecies, double, double *),
                    double (*dx4dt)(FourSpecies, double, double *),
                    double *dxdt_parameters, double dt, int steps) {
    /* Uses Euler's method to approximate the solution of an arbitrary four species ODE
                    dx1/dt = f1(x1, x2, x3, x4, t, parameters),
                    ...
                    dx4/dt = f4(x1, x2, x3, x4, t, parameters)
    x: a pointer to a TwoSpecies array of size steps to be filled with x values solved 
    t: a pointer to an array of size steps to be filled with t values
    dx1dt: a pointer to the function f1(x1, x2, x3, x4, t, parameters)
    dxdt_parameters: an array that holds parameters for the functions f1, ...
    dt: time step
    steps: number of steps to iterate Euler's method */

    int i;
    for (i = 0; i < steps - 1; i++) {
        double f1 = (*dx1dt)(*(x+i), *(t+i), dxdt_parameters);
        double f2 = (*dx2dt)(*(x+i), *(t+i), dxdt_parameters);
        double f3 = (*dx3dt)(*(x+i), *(t+i), dxdt_parameters);
        double f4 = (*dx4dt)(*(x+i), *(t+i), dxdt_parameters);
        *(x+i+1) = euler_step_fourspecies(*(x+i), f1, f2, f3, f4, dt);
        *(t+i+1) = *(t+i) + dt;
    }
}

FourSpecies euler_step_fourspecies(FourSpecies x, double f1, double f2, double f3, double f4, double dt) {
    double new_x1 = x.x1 + f1 * dt;
    double new_x2 = x.x2 + f2 * dt;
    double new_x3 = x.x3 + f3 * dt;
    double new_x4 = x.x4 + f4 * dt;
    FourSpecies new_x = {new_x1, new_x2, new_x3, new_x4};
    return new_x;
}