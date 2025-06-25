#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>

#include <math.h>
#include <iostream>

#include <sundials/sundials_core.h>     // Provides core SUNDIALS types
#include <cvode/cvode.h>                // CVODE provides linear multistep methods
#include <nvector/nvector_serial.h>     // access to serial N_Vector
#include <sunlinsol/sunlinsol_dense.h>  // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h>  // access to dense SUNmatrix

//#include <cvode/cvode.h>
//#include <nvector/nvector_serial.h>
//#include <cvode/cvode_dense.h>

#include "de.h"


// SUNDIALS context
static SUNContext sunctx = NULL;

#define SCALE 12.5

#define ORIGIN         0
#define CENTER_OF_MASS 1

int connections[2*12] = {
    0, 1,
    0, 2,
    0, 3,
    0, 4,
    0, 5,
    0, 6,
    1, 2,
    2, 3,
    3, 4,
    4, 5,
    5, 6,
    6, 1
};



class Playback : public Fl_Gl_Window {

    int frame;
    double frame_period;

    int center;

    xparams_t params;
    N_Vector state;

    void *cvode_mem;
    sunrealtype tau;
    sunrealtype dtau;
    sunrealtype tau1;

    int crashed;

    //
    // Draw the system.
    //
    void draw()
    {
        double ox, oy;
        double xx, yy;
        int idx;

        if (!valid()) {
            valid(1);
            glLoadIdentity();
            glViewport(0, 0, w(), h());
        }
        glClear(GL_COLOR_BUFFER_BIT);

        ox = 0.0;
        oy = 0.0;
        if (center == CENTER_OF_MASS) {
            for (idx = 0; idx < params.num_points; ++idx) {
                ox += NV_Ith_S(state, 2*idx);
                oy += NV_Ith_S(state, 2*idx + 1);
            }
            ox /= params.num_points;
            oy /= params.num_points;
        }
        glPushMatrix();
            // Draw the connections. Color-code based on whether the spring is
            // stretched or compressed (relative to the natural length params.L).
            for (idx = 0; idx < params.num_connections; ++idx) {
                int i, j;
                double red, blue;
                double line_base_color = 0.6;
                i = params.connections[2*idx];
                j = params.connections[2*idx+1];

                double x1 = NV_Ith_S(state, 2*i);
                double y1 = NV_Ith_S(state, 2*i+1);
                double x2 = NV_Ith_S(state, 2*j);
                double y2 = NV_Ith_S(state, 2*j+1);
                double dist = hypot(x2 - x1, y2 - y1);
                if (dist <= params.L) {
                    red = line_base_color + (1 - line_base_color)*tanh(25*(params.L - dist)/params.L);
                    blue = line_base_color;
                }
                else if (dist > params.L) {
                    blue = line_base_color + (1 - line_base_color)*tanh(25*(dist - params.L)/params.L);
                    red = line_base_color;
                }
                glColor3f(red, line_base_color, blue);

                glBegin(GL_LINES);
                xx = (x1 - ox)/SCALE;
                yy = (y1 - oy)/SCALE;
                glVertex2f(xx, yy);
                xx = (x2 - ox)/SCALE;
                yy = (y2 - oy)/SCALE;
                glVertex2f(xx, yy);
                glEnd();
            }

            // Draw the point masses as white dots.
            glColor3f(1.0, 1.0, 1.0);
            glPointSize(3.0);
            glBegin(GL_POINTS);
            for (int i = 0; i < params.num_points; ++i) {
                double x = (NV_Ith_S(state, 2*i) - ox)/SCALE;
                double y = (NV_Ith_S(state, 2*i+1) - oy)/SCALE;
                glVertex2f(x, y);
            }
            glEnd();

            // Draw the disk in the center.
            glColor3f(0.75, 0.75, 0.0);
            glBegin(GL_POLYGON);
            for (int i = 0; i < 13; i++) {
                double theta = 2*M_PI * i / 13;
                xx = (params.r0*cos(theta) - ox)/SCALE;
                yy = (params.r0*sin(theta) - oy)/SCALE;
                glVertex2f(xx, yy);
            }
            glEnd();
        glPopMatrix();

        // Advance frame counter
        ++frame;
    }

    //
    // Called repeatedly to redraw the window.
    //
    static void Timer_CB(void *userdata)
    {
        int flag;
        Playback *pb = (Playback*)userdata;

        if (pb->cvode_mem == NULL) {
            std::cerr << "pb->cvode_mem is NULL!\n";
        }
        flag = CVode(pb->cvode_mem, pb->tau + pb->dtau, pb->state,
                     &(pb->tau), CV_NORMAL);
        if (flag == CV_ROOT_RETURN) {
            pb->crashed = 1;
        }
        pb->redraw();
        if (!pb->crashed) {
            Fl::repeat_timeout(pb->frame_period, Timer_CB, userdata);
        }
    }

public:

    int handle(int e)
    {
        char c;

        if (e == FL_FOCUS || e == FL_UNFOCUS) {
            return 1;
        }
        if (e == FL_KEYBOARD) {
            c = Fl::event_text()[0];
            std::cerr << Fl::event_key() << " '" << c << "'\n";
            if (c == '+') {
                dtau *= 2.0;
            }
            else if (c == '-') {
                dtau *= 0.5;
            }
            else if (c == 'o') {
                center = ORIGIN;
            }
            else if (c == 'c') {
                center = CENTER_OF_MASS;
            }
        }
        return(Fl_Gl_Window::handle(e));
    }

    // Constructor
    Playback(int X, int Y, int W, int H, const char*L=0) : Fl_Gl_Window(X,Y,W,H,L)
    {
        int retval;
        int flag;
        double *p;

        std::cerr << "Playback() start\n";

        /* Create the SUNDIALS context */
        retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
        if (retval) {
            std::cerr << "SUNContext_Create() failed.\n";
            return;
        }

        frame = 0;
        frame_period = 1.0/24.0;

        center = ORIGIN;

        params.k = 2.5;
        params.L = 1.5;
        params.b = 0.5;
        params.g = 8.0;
        params.r0 = 0.25;
        params.num_points = 7;
        params.num_connections = 12;
        params.connections = connections;

        state = N_VNew_Serial(4*params.num_points, sunctx);
        p = N_VGetArrayPointer(state);
        hex_ics(9.0, 9.0, params.L, p);

        double v0 = 0.45;
        NV_Ith_S(state, 14) =  v0;
        NV_Ith_S(state, 15) = -v0;
        NV_Ith_S(state, 16) =  v0;
        NV_Ith_S(state, 17) = -v0;
        NV_Ith_S(state, 18) =  v0;
        NV_Ith_S(state, 19) = -v0;
        NV_Ith_S(state, 20) =  v0;
        NV_Ith_S(state, 21) = -v0;
        NV_Ith_S(state, 22) =  v0;
        NV_Ith_S(state, 23) = -v0;
        NV_Ith_S(state, 24) =  v0;
        NV_Ith_S(state, 25) = -v0;
        NV_Ith_S(state, 26) = 0.1*v0;
        NV_Ith_S(state, 27) = -0.1*v0;

        tau = SUN_RCONST(0.0);
        dtau = SUN_RCONST(0.125);
        tau1 = SUN_RCONST(2500.0);

        cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
        if (cvode_mem == NULL) {
            std::cerr << "darn it\n";
        }

        flag = CVodeInit(cvode_mem, de, tau, state);
        flag = CVodeSStolerances(cvode_mem, 1e-10, 1e-12);
        flag = CVodeSetUserData(cvode_mem, &params);
        flag = CVodeSetMaxNumSteps(cvode_mem, 500000);

        // Create dense SUNMatrix for use in linear solves
        SUNMatrix A = SUNDenseMatrix(4*params.num_points, 4*params.num_points, sunctx);
        if (!A) {
            std::cerr << "SUNDenseMatrix() failed.\n";
        }

        // Create dense SUNLinearSolver object
        SUNLinearSolver LS = SUNLinSol_Dense(state, A, sunctx);
        if (!LS) {
            std::cerr << "SUNLinSol_Dense() failed.\n";
        }

        // Attach the matrix and linear solver to CVODE
        retval = CVodeSetLinearSolver(cvode_mem, LS, A);
        if (retval) {
            std::cerr << "CVodeSetLinearSolver() failed.\n";
        }

        flag = CVodeSetStopTime(cvode_mem, tau1);
        flag = CVodeRootInit(cvode_mem, params.num_points, collision);

        crashed = 0;

        Fl::add_timeout(frame_period, Timer_CB, (void*)this);
        end();

        std::cerr << "Playback() returning\n";
    }
};


int main()
{
    Fl_Window win(720, 720);
    Playback playback(10, 10, win.w()-20, win.h()-20);
    win.resizable(&playback);
    win.show();
    return(Fl::run());
}
