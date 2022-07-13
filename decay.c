/* Ravi Chepuri */

#include <stdio.h>
#include <stdlib.h> // includes the strtof function to convert strings to floats
#include "euler_solver.h"
#include "rk2_solver.h"
#include "rk4_solver.h"


#define INIT_COND 100
#define DT0 0.5
#define DT1 0.2
#define DT2 0.05
#define MAX_TIME 6.0

double f(double x, double t, double *parameters);


int main(int argc, char *argv[]) {
    /* Given a list of time-steps provided as arguments to the executable at the command line, calculate single-species decay with each time-step. */

    double fparams[] = {1.}; //time constant
    float dt[argc-1];
    int steps[argc-1];
    int i, j;
    for (j = 0; j < argc-1; j++) {

        dt[j] = strtof(argv[j+1], NULL);
        steps[j] = (int) (MAX_TIME / dt[j]);

        double n_euler[steps[j]], n_rk2[steps[j]], n_rk4[steps[j]], t[steps[j]];
        n_euler[0] = INIT_COND, n_rk2[0] = INIT_COND, n_rk4[0] = INIT_COND;

        t[0] = 0.;
        euler_solver(n_euler, t, &f, fparams, dt[j], steps[j]);
        t[0] = 0.;
        rk2_solver(n_rk2, t, &f, fparams, dt[j], steps[j]);
        t[0] = 0.;
        rk4_solver(n_rk4, t, &f, fparams, dt[j], steps[j]);

        char filename[40];
        snprintf(filename, 40, "./decay_data/%s.dat", argv[j+1]); // store the appropriate filename in the string "filename"
        FILE* fp = fopen(filename, "w+");
        fprintf(fp, "t,euler,rk2,rk4,\n");
        for (i = 0; i < steps[j]; i++) {
            fprintf(fp, "%f,%f,%f,%f,\n", t[i], n_euler[i], n_rk2[i], n_rk4[i]);
        }
        fclose(fp);
    }

    // /* Euler method */

    // double n0_euler[steps0], t0_euler[steps0]; // for time step DT0
    // double n1_euler[steps1], t1_euler[steps1]; // for time step DT1
    // double n2_euler[steps2], t2_euler[steps2]; // for time step DT2
    // n0_euler[0] = 100.0, t0_euler[0] = 0.0;
    // n1_euler[0] = 100.0, t1_euler[0] = 0.0;
    // n2_euler[0] = 100.0, t2_euler[0] = 0.0;
    // euler_solver(n0_euler, t0_euler, &f, fparams, DT0, steps0);
    // euler_solver(n1_euler, t1_euler, &f, fparams, DT1, steps1);
    // euler_solver(n2_euler, t2_euler, &f, fparams, DT2, steps2);
    // FILE* fp0_euler = fopen("./data/euler_dt0.dat", "w+");
    // FILE* fp1_euler = fopen("./data/euler_dt1.dat", "w+");
    // FILE* fp2_euler = fopen("./data/euler_dt2.dat", "w+");
    // fprintf(fp0_euler, "t, N\n");
    // fprintf(fp1_euler, "t, N\n");
    // fprintf(fp2_euler, "t, N\n");
    // for (i = 0; i < steps0; i++) {
    //     fprintf(fp0_euler, "%f, %f\n", t0_euler[i], n0_euler[i]);
    // }
    // for (i = 0; i < steps1; i++) {
    //     fprintf(fp1_euler, "%f, %f\n", t1_euler[i], n1_euler[i]);
    // }
    // for (i = 0; i < steps2; i++) {
    //     fprintf(fp2_euler, "%f, %f\n", t2_euler[i], n2_euler[i]);
    // }
    // fclose(fp0_euler);
    // fclose(fp1_euler);
    // fclose(fp2_euler);
    
    // /* Runge-Kutta Order 2 method */

    // double n0_rk2[steps0], t0_rk2[steps0];
    // double n1_rk2[steps1], t1_rk2[steps1];
    // double n2_rk2[steps2], t2_rk2[steps2];
    // n0_rk2[0] = 100.0, t0_rk2[0] = 0.0;
    // n1_rk2[0] = 100.0, t1_rk2[0] = 0.0;
    // n2_rk2[0] = 100.0, t2_rk2[0] = 0.0;
    // rk2_solver(n0_rk2, t0_rk2, &f, fparams, DT0, steps0);
    // rk2_solver(n1_rk2, t1_rk2, &f, fparams, DT1, steps1);
    // rk2_solver(n2_rk2, t2_rk2, &f, fparams, DT2, steps2);
    // FILE* fp0_rk2 = fopen("./data/rk2_dt0.dat", "w+");
    // FILE* fp1_rk2 = fopen("./data/rk2_dt1.dat", "w+");
    // FILE* fp2_rk2 = fopen("./data/rk2_dt2.dat", "w+");
    // fprintf(fp0_rk2, "t, N\n");
    // fprintf(fp1_rk2, "t, N\n");
    // fprintf(fp2_rk2, "t, N\n");
    // for (i = 0; i < steps0; i++) {
    //     fprintf(fp0_rk2, "%f, %f\n", t0_rk2[i], n0_rk2[i]);
    // }
    // for (i = 0; i < steps1; i++) {
    //     fprintf(fp1_rk2, "%f, %f\n", t1_rk2[i], n1_rk2[i]);
    // }
    // for (i = 0; i < steps2; i++) {
    //     fprintf(fp2_rk2, "%f, %f\n", t2_rk2[i], n2_rk2[i]);
    // }
    // fclose(fp0_rk2);
    // fclose(fp1_rk2);
    // fclose(fp2_rk2);

    // /* Runge Kutta Order 4 method */

    // double n0_rk4[steps0], t0_rk4[steps0];
    // double n1_rk4[steps1], t1_rk4[steps1];
    // double n2_rk4[steps2], t2_rk4[steps2];
    // n0_rk4[0] = 100.0, t0_rk4[0] = 0.0;
    // n1_rk4[0] = 100.0, t1_rk4[0] = 0.0;
    // n2_rk4[0] = 100.0, t2_rk4[0] = 0.0;
    // rk4_solver(n0_rk4, t0_rk4, &f, fparams, DT0, steps0);
    // rk4_solver(n1_rk4, t1_rk4, &f, fparams, DT1, steps1);
    // rk4_solver(n2_rk4, t2_rk4, &f, fparams, DT2, steps2);
    // FILE* fp0_rk4 = fopen("./data/rk4_dt0.dat", "w+");
    // FILE* fp1_rk4 = fopen("./data/rk4_dt1.dat", "w+");
    // FILE* fp2_rk4 = fopen("./data/rk4_dt2.dat", "w+");
    // fprintf(fp0_rk4, "t, N\n");
    // fprintf(fp1_rk4, "t, N\n");
    // fprintf(fp2_rk4, "t, N\n");
    // for (i = 0; i < steps0; i++) {
    //     fprintf(fp0_rk4, "%f, %f\n", t0_rk4[i], n0_rk4[i]);
    // }
    // for (i = 0; i < steps1; i++) {
    //     fprintf(fp1_rk4, "%f, %f\n", t1_rk4[i], n1_rk4[i]);
    // }
    // for (i = 0; i < steps2; i++) {
    //     fprintf(fp2_rk4, "%f, %f\n", t2_rk4[i], n2_rk4[i]);
    // }
    // fclose(fp0_rk4);
    // fclose(fp1_rk4);
    // fclose(fp2_rk4);
}


double f(double x, double t, double *parameters) {
    return -x / parameters[0];
}