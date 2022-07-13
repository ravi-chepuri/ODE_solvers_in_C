/* Ravi Chepuri */

#include "two_species.h"
#include "four_species.h"

void rk2_solver(double *x, double *t, double (*dxdt)(double, double, double *),
                    double* dxdt_parameters, double dt, int steps);

void rk2_solver_twospecies(TwoSpecies *x, double *t,
                    double (*dxAdt)(TwoSpecies, double, double *),
                    double (*dxBdt)(TwoSpecies, double, double *),
                    double *dxdt_parameters, double dt, int steps);

void rk2_solver_fourspecies(FourSpecies *x, double *t,
                    double (*dx1dt)(FourSpecies, double, double *),
                    double (*dx2dt)(FourSpecies, double, double *),
                    double (*dx3dt)(FourSpecies, double, double *),
                    double (*dx4dt)(FourSpecies, double, double *),
                    double *dxdt_parameters, double dt, int steps);