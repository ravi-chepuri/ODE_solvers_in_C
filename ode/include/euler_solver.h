/* Ravi Chepuri */

#include "two_species.h"
#include "four_species.h"

double euler_step(double x, double f, double dt);

void euler_solver(double *x, double *t, double (*dxdt)(double, double, double *),
                    double* dxdt_parameters, double dt, int steps);


void euler_solver_twospecies(TwoSpecies *, double *,
                    double (*dxAdt)(TwoSpecies, double, double *),
                    double (*dxBdt)(TwoSpecies, double, double *),
                    double *, double, int);

TwoSpecies euler_step_twospecies(TwoSpecies, double, double, double);


void euler_solver_fourspecies(FourSpecies *, double *,
                    double (*dx1dt)(FourSpecies, double, double *),
                    double (*dx2dt)(FourSpecies, double, double *),
                    double (*dx3dt)(FourSpecies, double, double *),
                    double (*dx4dt)(FourSpecies, double, double *),
                    double *, double, int);

FourSpecies euler_step_fourspecies(FourSpecies, double, double, double, double, double);