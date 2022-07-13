/* Ravi Chepuri */

#include <stdio.h>
#include "rk4_solver.h"

void rk4_solver(double *x, double *t, double (*dxdt)(double, double, double *),
                    double* dxdt_parameters, double dt, int steps) {
    /* Uses a second order Runge-Kutta method to approximate the solution of an arbitrary ODE
                        dx/dt = f(x, t, parameters)
        x: a pointer to an array of size steps to be filled with x values solved 
        t: a pointer to an array of size steps to be filled with t values
        dxdt: a pointer to the function f(x, t, parameters)
        dt: time step
        steps: number of steps to iterate */

    int i;
    double t1p, t2p, t3p, t4p; // (t_1)', etc. in the O(4) Runge-Kutta method
    double x1p, x2p, x3p, x4p; // (x_1)', etc.
    for (i = 0; i < steps - 1; i++) {
        /* Equations on slide 16 of Lecture 5 */
        x1p = *(x+i);
        t1p = *(t+i);
        x2p = *(x+i) + (*dxdt)(x1p, t1p, dxdt_parameters) * dt / 2;
        t2p = *(t+i) + dt / 2;
        x3p = *(x+i) + (*dxdt)(x2p, t2p, dxdt_parameters) * dt / 2;
        t2p = *(t+i) + dt / 2;
        x4p = *(x+i) + (*dxdt)(x3p, t3p, dxdt_parameters) * dt;
        t4p = *(t+i) + dt;
        *(x+i+1) = *(x+i) + (dt/6.) * (   1 * (*dxdt)(x1p, t1p, dxdt_parameters)
                                        + 2 * (*dxdt)(x2p, t2p, dxdt_parameters)
                                        + 2 * (*dxdt)(x3p, t3p, dxdt_parameters)
                                        + 1 * (*dxdt)(x4p, t4p, dxdt_parameters) );
        *(t+i+1) = *(t+i) + dt;
    }
}


void rk4_solver_twospecies(TwoSpecies *x, double *t,
                    double (*dxAdt)(TwoSpecies, double, double *),
                    double (*dxBdt)(TwoSpecies, double, double *),
                    double *dxdt_parameters, double dt, int steps) {
    /* Uses a fourth order Runge Kutta method to approximate the solution of an arbitrary two species ODE
                        dx_A/dt = f_A(x_A, x_B, t, parameters)
                        dx_B/dt = f_B(x_A, x_B, t, parameters)
        x: a pointer to a TwoSpecies array of size steps to be filled with x values solved 
        t: a pointer to an array of size steps to be filled with t values
        dxAdt: a pointer to the function f_A(x_A, x_B, t, parameters)
        dxBdt: a pointer to the function f_B(x_A, x_B, t, parameters)
        dxdt_parameters: an array that holds parameters for the function f(x, t, parameters)
        dt: time step
        steps: number of steps to iterate Euler's method */
    
    double x1pA, x1pB, x2pA, x2pB, x3pA, x3pB, x4pA, x4pB;
    TwoSpecies x1p, x2p, x3p, x4p;
    double t1p, t2p, t3p, t4p;
    double f1A, f1B, f2A, f2B, f3A, f3B, f4A, f4B;
    double new_xA, new_xB;
    TwoSpecies new_x;
    int i;
    
    for (i = 0; i < steps - 1; i++) {

        /* Equations on slide 16 of Lecture 5 */

        x1pA = (x+i)->Na;
        x1pB = (x+i)->Nb;
        x1p.Na = x1pA;
        x1p.Nb = x1pB;
        t1p = *(t+i);

        f1A = (*dxAdt)(x1p, t1p, dxdt_parameters);
        f1B = (*dxBdt)(x1p, t1p, dxdt_parameters);

        x2pA = (x+i)->Na + 0.5 * f1A * dt;
        x2pB = (x+i)->Nb + 0.5 * f1B * dt;
        x2p.Na = x2pA;
        x2p.Nb = x2pB;
        t2p = *(t+i) + 0.5 * dt;

        f2A = (*dxAdt)(x2p, t2p, dxdt_parameters);
        f2B = (*dxBdt)(x2p, t2p, dxdt_parameters);

        x3pA = (x+i)->Na + 0.5 * f2A * dt;
        x3pB = (x+i)->Nb + 0.5 * f2B * dt;
        x3p.Na = x3pA;
        x3p.Nb = x3pB;
        t3p = *(t+i) + 0.5 * dt;

        f3A = (*dxAdt)(x3p, t3p, dxdt_parameters);
        f3B = (*dxBdt)(x3p, t3p, dxdt_parameters);

        x4pA = (x+i)->Na + f3A * dt;
        x4pB = (x+i)->Nb + f3B * dt;
        x4p.Na = x4pA;
        x4p.Nb = x4pB;
        t4p = *(t+i) + dt;

        f4A = (*dxAdt)(x4p, t4p, dxdt_parameters);
        f4B = (*dxBdt)(x4p, t4p, dxdt_parameters);

        new_xA = (x+i)->Na + 1./6. * (f1A + 2*f2A + 2*f3A + f4A) * dt;
        new_xB = (x+i)->Nb + 1./6. * (f1B + 2*f2B + 2*f3B + f4B) * dt;
        new_x.Na = new_xA;
        new_x.Nb = new_xB;

        *(x+i+1) = new_x;
        *(t+i+1) = *(t+i) + dt;
    }
}


void rk4_solver_fourspecies(FourSpecies *x, double *t,
                    double (*dx1dt)(FourSpecies, double, double *),
                    double (*dx2dt)(FourSpecies, double, double *),
                    double (*dx3dt)(FourSpecies, double, double *),
                    double (*dx4dt)(FourSpecies, double, double *),
                    double *dxdt_parameters, double dt, int steps) {
    /* Uses a fourth order Runge Kutta method to approximate the solution of an arbitrary four species ODE
                dx1/dt = f1(x1, x2, x3, x4, t, parameters),
                ...
                dx4/dt = f4(x1, x2, x3, x4, t, parameters),
    x: a pointer to a TwoSpecies array of size steps to be filled with x values solved 
    t: a pointer to an array of size steps to be filled with t values
    dx1dt: a pointer to the function f1(x1, x2, x3, x4, t, parameters)
    dxdt_parameters: an array that holds parameters for the functions f1, ...
    dt: time step
    steps: number of steps to iterate */

    double  x1p1, x1p2, x1p3, x1p4, // 1st component of x_1', 2nd component of x_1', etc.
            x2p1, x2p2, x2p3, x2p4, // 1st component of x_2', etc.
            x3p1, x3p2, x3p3, x3p4, 
            x4p1, x4p2, x4p3, x4p4;
    FourSpecies x1p, x2p, x3p, x4p; // x_1', etc.
    double t1p, t2p, t3p, t4p; //t_1', etc.
    double  f11, f12, f13, f14, // 1st component of f(x_1', t_1'), 2nd component of ..., etc.
            f21, f22, f23, f24, // 1st component of f(x_2', t_2'), etc.
            f31, f32, f33, f34,
            f41, f42, f43, f44;
    double new_x1, new_x2, new_x3, new_x4;
    FourSpecies new_x;
    int i;

    for (i = 0; i < steps - 1; i++) {

        /* Equations on slide 16 of Lecture 5 */

        x1p1 = (x+i)->x1;
        x1p2 = (x+i)->x2;
        x1p3 = (x+i)->x3;
        x1p4 = (x+i)->x4;
        x1p.x1 = x1p1;
        x1p.x2 = x1p2;
        x1p.x3 = x1p3;
        x1p.x4 = x1p4;
        t1p = *(t+i);

        f11 = (*dx1dt)(x1p, t1p, dxdt_parameters);
        f12 = (*dx2dt)(x1p, t1p, dxdt_parameters);
        f13 = (*dx3dt)(x1p, t1p, dxdt_parameters);
        f14 = (*dx4dt)(x1p, t1p, dxdt_parameters);

        x2p1 = (x+i)->x1 + 0.5 * f11 * dt;
        x2p2 = (x+i)->x2 + 0.5 * f12 * dt;
        x2p3 = (x+i)->x3 + 0.5 * f13 * dt;
        x2p4 = (x+i)->x4 + 0.5 * f14 * dt;
        x2p.x1 = x2p1;
        x2p.x2 = x2p2;
        x2p.x3 = x2p3;
        x2p.x4 = x2p4;
        t2p = *(t+i) + 0.5 * dt;

        f21 = (*dx1dt)(x2p, t2p, dxdt_parameters);
        f22 = (*dx2dt)(x2p, t2p, dxdt_parameters);
        f23 = (*dx3dt)(x2p, t2p, dxdt_parameters);
        f24 = (*dx4dt)(x2p, t2p, dxdt_parameters);

        x3p1 = (x+i)->x1 + 0.5 * f21 * dt;
        x3p2 = (x+i)->x2 + 0.5 * f22 * dt;
        x3p3 = (x+i)->x3 + 0.5 * f23 * dt;
        x3p4 = (x+i)->x4 + 0.5 * f24 * dt;
        x3p.x1 = x3p1;
        x3p.x2 = x3p2;
        x3p.x3 = x3p3;
        x3p.x4 = x3p4;
        t3p = *(t+i) + 0.5 * dt;

        f31 = (*dx1dt)(x3p, t3p, dxdt_parameters);
        f32 = (*dx2dt)(x3p, t3p, dxdt_parameters);
        f33 = (*dx3dt)(x3p, t3p, dxdt_parameters);
        f34 = (*dx4dt)(x3p, t3p, dxdt_parameters);

        x4p1 = (x+i)->x1 + f31 * dt;
        x4p2 = (x+i)->x2 + f32 * dt;
        x4p3 = (x+i)->x3 + f33 * dt;
        x4p4 = (x+i)->x4 + f34 * dt;
        x4p.x1 = x4p1;
        x4p.x2 = x4p2;
        x4p.x3 = x4p3;
        x4p.x4 = x4p4;
        t4p = *(t+i) + dt;

        f41 = (*dx1dt)(x4p, t4p, dxdt_parameters);
        f42 = (*dx2dt)(x4p, t4p, dxdt_parameters);
        f43 = (*dx3dt)(x4p, t4p, dxdt_parameters);
        f44 = (*dx4dt)(x4p, t4p, dxdt_parameters);

        new_x1 = (x+i)->x1 + 1./6. * (f11 + 2*f21 + 2*f31 + f41) * dt;
        new_x2 = (x+i)->x2 + 1./6. * (f12 + 2*f22 + 2*f32 + f42) * dt;
        new_x3 = (x+i)->x3 + 1./6. * (f13 + 2*f23 + 2*f33 + f43) * dt;
        new_x4 = (x+i)->x4 + 1./6. * (f14 + 2*f24 + 2*f34 + f44) * dt;

        new_x.x1 = new_x1;
        new_x.x2 = new_x2;
        new_x.x3 = new_x3;
        new_x.x4 = new_x4;

        *(x+i+1) = new_x;
        *(t+i+1) = *(t+i) + dt;
    }
}   