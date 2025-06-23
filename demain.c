#include <stdio.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>

#include "de.h"

int main(int argc, char *argv[])
{
    int flag;
    int j;

    params_t p;

    p.k = 0.25;
    p.L = 1.0;
    p.b = 4.0;
    p.g = 5.0;

    /* Initial conditions */
    N_Vector w;
    w = N_VNew_Serial(12);
    X(w,0) = 6.0;
    Y(w,0)= 8.2;
    X(w,1) = 8.2;
    Y(w,1) = 8.0;
    X(w,2) = 8.2;
    Y(w,2) = 8.2;
    U(w,0) = 0.4;
    V(w,0) = -0.4;
    U(w,1) = 0.1;
    V(w,1) = -0.25;
    U(w,2) = 0.5;
    V(w,2) = -0.25;

    /* For non-stiff problems:   */
    void *cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
    /* For stiff problems:       */
    //void *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    realtype t = RCONST(0.0);
    flag = CVodeInit(cvode_mem, de3, t, w);
    flag = CVodeSStolerances(cvode_mem, 1e-10, 1e-12);
    flag = CVodeSetUserData(cvode_mem, &p);
    flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
    flag = CVDense(cvode_mem, 12);

    realtype dt = RCONST(0.25);
    realtype t1 = RCONST(2500.0);
    flag = CVodeSetStopTime(cvode_mem, t1);

    /* Print the solution at the current time */
    printf("%.8e", t);
    for (j = 0; j < 12; ++j) {
        printf(" %.8e", NV_Ith_S(w, j));
    }
    printf("\n");

    while (t < t1) {
        /* Advance the solution */
        //flag = CVode(cvode_mem, t+dt, w, &t, CV_ONE_STEP);
        flag = CVode(cvode_mem, t+dt, w, &t, CV_NORMAL);
        if (flag != CV_SUCCESS && flag != CV_TSTOP_RETURN) {
            fprintf(stderr, "flag=%d\n", flag);
            break;
        }
        /* Print the solution at the current time */
        printf("%.8e", t);
        for (j = 0; j < 12; ++j) {
            printf(" %.8e", NV_Ith_S(w, j));
        }
        printf("\n");
        //t = t + dt;
    }
    N_VDestroy_Serial(w);
    CVodeFree(&cvode_mem);
}
