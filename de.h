#ifndef _DE_H_
#define _DE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <nvector/nvector_serial.h> // access to serial N_Vector

typedef struct _params {
    double k, L, b, g;
} params_t;

typedef struct _xparams {
    double k, L, b, g;
    double r0;
    /* num_points is the number of points */
	int num_points;
	/* num_connections is the number of connections between points */
    int num_connections;
    /* connections is an array of length 2*n */
    int *connections;
} xparams_t;

typedef struct _rigid_hex_params {
	double L, g;
	double r0;
} rigid_hex_params_t;

#define X(w,i)  NV_Ith_S(w, 2*(i))
#define Y(w,i)  NV_Ith_S(w, 2*(i)+1)
#define U(w,i)  NV_Ith_S(w, 2*(i)+6)
#define V(w,i)  NV_Ith_S(w, 2*(i)+7)

int de3(sunrealtype t, N_Vector w, N_Vector f, void *params);
int de(sunrealtype t, N_Vector w, N_Vector f, void *params);
int collision(sunrealtype t, N_Vector y, sunrealtype *gout, void *user_data);
int hex_ics(double cx, double cy, double L, double *p);
int de_rigid_hex(sunrealtype t, N_Vector w, N_Vector f, void *params);

#ifdef __cplusplus
}
#endif

#endif
