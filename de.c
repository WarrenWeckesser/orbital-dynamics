#include <sundials/sundials_core.h> // Provides core SUNDIALS types

#include <math.h>
#include <stdio.h>
#include "de.h"

#ifdef __cplusplus
extern "C" {
#endif

//
// Force of the spring, given separation r, spring constant k,
// and natural length L.
//
double spring_force(double r, double k, double L)
{
    double rho = L/r;
    return -k*r*(1 - rho)*(1 + rho + rho*rho);
    // return k*r*(1 - rho);
    // return k*(L - r);
}


//
//  The vector field for three point masses connected by springs.
//
// (This function has not been maintained lately.)
//
//  w holds [x0, y0, x1, y1, x2, y2,
//            dx0/dt, dy0/dt, dx1/dt, dy1/dt, dx2/dt, dy2/dt]
//

int de3(sunrealtype t, N_Vector w, N_Vector f, void *params)
{
    params_t *p = params;
    int indices[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    int idx;

    for (idx = 0; idx < 3; ++idx) {
        X(f,idx) = U(w,idx);
        Y(f,idx) = V(w,idx);
    }

    for (idx = 0; idx < 3; ++idx) {
        int i, j;
        double uvec[2], dist, fvec[2];
        double relvel[2], s;
        double force;

        i = indices[idx][0];
        j = indices[idx][1];

        uvec[0] = X(w,j) - X(w,i);
        uvec[1] = Y(w,j) - Y(w,i);
        dist = hypot(uvec[0], uvec[1]);
        uvec[0] /= dist;
        uvec[1] /= dist;

        force = spring_force(dist, p->k, p->L);
        fvec[0] = force*uvec[0];
        fvec[1] = force*uvec[1];
        U(f,i) += fvec[0];
        V(f,i) += fvec[1];
        U(f,j) -= fvec[0];
        V(f,j) -= fvec[1];

        relvel[0] = U(w,j) - U(w,i);
        relvel[1] = V(w,j) - V(w,i);
        s = relvel[0]*uvec[0] + relvel[1]*uvec[1];
        U(f,i) += p->b * s * uvec[0];
        V(f,i) += p->b * s * uvec[1];
        U(f,j) -= p->b * s * uvec[0];
        V(f,j) -= p->b * s * uvec[1];
    }

    if (p->g > 0) {
        for (idx = 0; idx < 3; ++idx) {
            double r = hypot(X(w, idx), Y(w, idx));
            if (r < 0.0001) {
                r = 0.0001;
            }
            double r3 = r*r*r;
            U(f,idx) += -p->g * X(w,idx) / r3;
            V(f,idx) += -p->g * Y(w,idx) / r3;
        }
    }

    return 0;
}


//
// The vector field for a system of point masses in a plane connected
// by springs.
//
//
//  w holds [x0, y0, x1, y1, x2, y2, ...,
//           dx0/dt, dy0/dt, dx1/dt, dy1/dt, dx2/dt, dy2/dt, ...]
//


int de(sunrealtype t, N_Vector w, N_Vector f, void *params)
{
    xparams_t *p = params;
    int num_points;
    int idx;
    int num_connections;
    int *connections;

    num_points = p->num_points;
    num_connections = p->num_connections;
    connections = p->connections;

    for (idx = 0; idx < num_points; ++idx) {
        NV_Ith_S(f, 2*idx) = NV_Ith_S(w, 2*num_points + 2*idx);
        NV_Ith_S(f, 2*idx+1) = NV_Ith_S(w, 2*num_points + 2*idx + 1);
        NV_Ith_S(f, 2*num_points + 2*idx) = 0.0;
        NV_Ith_S(f, 2*num_points + 2*idx + 1) = 0.0;
    }

    for (idx = 0; idx < num_connections; ++idx) {
        int i, j, ii, jj;
        double uvec[2], dist, fvec[2];
        double relvel[2], s;
        double force;

        i = connections[2*idx];
        j = connections[2*idx + 1];

        // uvec is a unit vector pointing from point i to point j.
        uvec[0] = NV_Ith_S(w, 2*j) - NV_Ith_S(w, 2*i);
        uvec[1] = NV_Ith_S(w, 2*j+1) - NV_Ith_S(w, 2*i+1);
        dist = hypot(uvec[0], uvec[1]);
        uvec[0] /= dist;
        uvec[1] /= dist;

        //
        // Compute the contribution of the spring forces to the system of equations.
        //
        force = spring_force(dist, p->k, p->L);
        fvec[0] = force*uvec[0];
        fvec[1] = force*uvec[1];
        ii = 2*num_points + 2*i;
        jj = 2*num_points + 2*j;
        NV_Ith_S(f, ii)     -= fvec[0];
        NV_Ith_S(f, ii + 1) -= fvec[1];
        NV_Ith_S(f, jj)     += fvec[0];
        NV_Ith_S(f, jj + 1) += fvec[1];

        //
        // Compute the contribution of the friction forces to the system of equations.
        //
        // relvel is the velocity of point j relative to point i.
        relvel[0] = NV_Ith_S(w, jj) - NV_Ith_S(w, ii);
        relvel[1] = NV_Ith_S(w, jj+1) - NV_Ith_S(w, ii+1);
        // s is the rate of change of the distance between points i and j.
        s = relvel[0]*uvec[0] + relvel[1]*uvec[1];
        // Friction force between points i and j.
        force = p->b * s;
        fvec[0] = force * uvec[0];
        fvec[1] = force * uvec[1];

        NV_Ith_S(f, ii)   += fvec[0];
        NV_Ith_S(f, ii+1) += fvec[1];
        NV_Ith_S(f, jj)   -= fvec[0];
        NV_Ith_S(f, jj+1) -= fvec[1];
    }

    if (p->g > 0) {
        // Compute the contribution of the gravitational forces to
        // the system of equations.
        for (idx = 0; idx < num_points; ++idx) {
            double xi, yi;
            xi = NV_Ith_S(w, 2*idx);
            yi = NV_Ith_S(w, 2*idx+1);
            double r = hypot(xi, yi);
            double r3 = r*r*r;
            NV_Ith_S(f, 2*num_points + 2*idx) += -p->g * xi / r3;
            NV_Ith_S(f, 2*num_points + 2*idx + 1) += -p->g * yi / r3;
        }
    }

    return 0;
}

int collision(sunrealtype t, N_Vector w, sunrealtype *gout, void *user_data)
{
    xparams_t *p = user_data;
    int num_points;
    int idx;
    double r0;

    num_points = p->num_points;
    r0 = p->r0;

    if (p->g > 0) {
        for (idx = 0; idx < num_points; ++idx) {
            double r = sqrt(NV_Ith_S(w, 2*idx) * NV_Ith_S(w, 2*idx) +
                            NV_Ith_S(w, 2*idx+1) * NV_Ith_S(w, 2*idx+1));
            gout[idx] = r - r0;
        }
    }
    else {
        for (idx = 0; idx < num_points; ++idx) {
            gout[idx] = 1.0;
        }   
    }
    return 0;
}


int hex_ics(double cx, double cy, double L, double *p)
{
    int i;

    p[0] = cx;
    p[1] = cy;
    for (i = 1; i < 7; ++i) {
        double theta = (i - 1) * 2*M_PI / 6.0;
        p[2*i] = cx + L*cos(theta);
        p[2*i+1] = cy + L*sin(theta);
    }
    return 0;
}


//
//  7 point masses arrange in a hexagonal pattern,
//  rigidly tied together.  State variables are (xc, yc, theta),
//  where (xc, yc) is the center of the hexagon.
//
int de_rigid_hex(sunrealtype t, N_Vector w, N_Vector f, void *params)
{
    double xc, yc, theta;
    double u, v, omega;
    int idx;
    double Ic;

    rigid_hex_params_t *p = params;

    // Moment of inertia, 6*m*r**2, with m = 1, and r is p->L.
    Ic = 6.0*(p->L)*(p->L);

    // Current position and angle of the 7 point hexagon.
    xc = NV_Ith_S(w, 0);
    yc = NV_Ith_S(w, 1);
    theta = NV_Ith_S(w, 2);
    // Current translational and angular velocities.
    u = NV_Ith_S(w, 3);
    v = NV_Ith_S(w, 4);
    omega = NV_Ith_S(w, 5);

    NV_Ith_S(f, 0) = u;     // dx/dt = u
    NV_Ith_S(f, 1) = v;     // dy/dt = v
    NV_Ith_S(f, 2) = omega; // dtheta/dt = omega

    NV_Ith_S(f, 3) = 0.0;
    NV_Ith_S(f, 4) = 0.0;
    NV_Ith_S(f, 5) = 0.0;
    for (idx = 0; idx < 7; ++idx) {
        double dx, dy, xi, yi, ri, ri3, fx, fy;

        if (idx < 6) {
            dx = (p->L)*cos(theta);
            dy = (p->L)*sin(theta);
        }
        else {
            // idx == 6 is the center point.
            dx = 0.0;
            dy = 0.0;
        }
        xi = xc + dx;
        yi = yc + dy;
        ri = sqrt(xi*xi + yi*yi);
        ri3 = ri*ri*ri;

        // Components of the gravitational force.
        fx = -p->g * xi / ri3;
        fy = -p->g * yi / ri3;

        NV_Ith_S(f, 3) += fx;
        NV_Ith_S(f, 4) += fy;
        NV_Ith_S(f, 5) += -dx*fy + dy*fx;  // torque about (xc,yc)
    }
    // Total mass is 7.
    NV_Ith_S(f, 3) /= 7.0;
    NV_Ith_S(f, 4) /= 7.0;
    NV_Ith_S(f, 5) /= Ic;

    return 0;
}

#ifdef __cplusplus
}
#endif
