/* Ravi Chepuri */

#include <stdio.h>
#include "rk2_solver.h"

void rk2_solver(double *x, double *t, double (*dxdt)(double, double, double *),
                    double* dxdt_parameters, double dt, int steps) {
    /* Uses a second order Runge-Kutta method to approximate the solution of an arbitrary ODE
                        dx/dt = f(x, t, parameters)
        x: a pointer to an array of size steps to be filled with x values solved 
        t: a pointer to an array of size steps to be filled with t values
        dxdt: a pointer to the function f(x, t, parameters)
        dt: time step
        steps: number of steps to iterate */

    int i;
    double tp; // t' in the O(2) Runge-Kutta method
    double xp; // x'
    for (i = 0; i < steps - 1; i++) {
        tp = *(t+i) + dt/2;
        xp = *(x+i) + (*dxdt)(*(x+i), *(t+i), dxdt_parameters) * dt/2;
        *(x+i+1) = *(x+i) + (*dxdt)(xp, tp, dxdt_parameters) * dt;
        *(t+i+1) = *(t+i) + dt;
    }
}

void rk2_solver_twospecies(TwoSpecies *x, double *t,
                    double (*dxAdt)(TwoSpecies, double, double *),
                    double (*dxBdt)(TwoSpecies, double, double *),
                    double *dxdt_parameters, double dt, int steps) {
    /* Uses a second order Runge Kutta method to approximate the solution of an arbitrary two species ODE
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
    double tp; // t' in the O(2) Runge-Kutta method
    double xAp; // x_A'
    double xBp; // x_B'
    TwoSpecies xp; // x' (a vector (x_A', x_B'))
    double new_xA;
    double new_xB;
    TwoSpecies new_x;

    for (i = 0; i < steps - 1; i++) {
        tp = *(t+i) + dt/2;
        xAp = (x+i)->Na + (*dxAdt)(*(x+i), *(t+i), dxdt_parameters) * dt/2;
        xBp = (x+i)->Nb + (*dxBdt)(*(x+i), *(t+i), dxdt_parameters) * dt/2;
        xp.Na = xAp;
        xp.Nb = xBp;
        new_xA = (x+i)->Na + (*dxAdt)(xp, tp, dxdt_parameters) * dt;
        new_xB = (x+i)->Nb + (*dxBdt)(xp, tp, dxdt_parameters) * dt;
        new_x.Na = new_xA;
        new_x.Nb = new_xB;
        *(x+i+1) = new_x;
        *(t+i+1) = *(t+i) + dt;
    }
}


void rk2_solver_fourspecies(FourSpecies *x, double *t,
                    double (*dx1dt)(FourSpecies, double, double *),
                    double (*dx2dt)(FourSpecies, double, double *),
                    double (*dx3dt)(FourSpecies, double, double *),
                    double (*dx4dt)(FourSpecies, double, double *),
                    double *dxdt_parameters, double dt, int steps) {
    /* Uses a second order Runge Kutta method to approximate the solution of an arbitrary four species ODE
                dx1/dt = f1(x1, x2, x3, x4, t, parameters),
                ...
                dx4/dt = f4(x1, x2, x3, x4, t, parameters),
    x: a pointer to a TwoSpecies array of size steps to be filled with x values solved 
    t: a pointer to an array of size steps to be filled with t values
    dx1dt: a pointer to the function f1(x1, x2, x3, x4, t, parameters)
    dxdt_parameters: an array that holds parameters for the functions f1, ...
    dt: time step
    steps: number of steps to iterate */

    int i;
    double f1, f2, f3, f4; // vector to store f(x, t)
    double x1p, x2p, x3p, x4p; // for x' vector
    FourSpecies xp; // x' vector
    double tp; // t'
    double f1p, f2p, f3p, f4p;
    double new_x1, new_x2, new_x3, new_x4;
    FourSpecies new_x; // new x value

    for (i = 0; i < steps - 1; i++) {
        f1 = (*dx1dt)(*(x+i), *(t+i), dxdt_parameters);
        f2 = (*dx2dt)(*(x+i), *(t+i), dxdt_parameters);
        f3 = (*dx3dt)(*(x+i), *(t+i), dxdt_parameters);
        f4 = (*dx4dt)(*(x+i), *(t+i), dxdt_parameters);
        tp = *(t+i) + dt/2;
        x1p = (x+i)->x1 + 0.5 * f1 * dt;
        x2p = (x+i)->x2 + 0.5 * f2 * dt;
        x3p = (x+i)->x3 + 0.5 * f3 * dt;
        x4p = (x+i)->x4 + 0.5 * f4 * dt;
        xp.x1 = x1p;
        xp.x2 = x2p;
        xp.x3 = x3p;
        xp.x4 = x4p;
        f1p = (*dx1dt)(xp, tp, dxdt_parameters);
        f2p = (*dx2dt)(xp, tp, dxdt_parameters);
        f3p = (*dx3dt)(xp, tp, dxdt_parameters);
        f4p = (*dx4dt)(xp, tp, dxdt_parameters);
        new_x1 = (x+i)->x1 + f1p * dt;
        new_x2 = (x+i)->x2 + f2p * dt;
        new_x3 = (x+i)->x3 + f3p * dt;
        new_x4 = (x+i)->x4 + f4p * dt;
        new_x.x1 = new_x1;
        new_x.x2 = new_x2;
        new_x.x3 = new_x3;
        new_x.x4 = new_x4;
        *(x+i+1) = new_x;
        *(t+i+1) = *(t+i) + dt;
    }
}