#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#include <math.h>
#include <iostream>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>

#include "de.h"

// Warning: Poorly designed C++ code ahead...

#define ORIGIN         0
#define CENTER_OF_MASS 1


class Playback : public Fl_Gl_Window {

    int frame;
    double frame_period;

    int center;

    rigid_hex_params_t params;
    N_Vector state;

    void *cvode_mem;
    realtype tau;
    realtype dtau;
    realtype tau1;

    int crashed;
    double scale;

    //
    // Draw the system.
    //
    void draw()
    {
        double ox, oy;
        double xc, yc, theta;
        double thetai, dx, dy, xx, yy;
        int idx;

        if (!valid()) {
            valid(1);
            glLoadIdentity();
            glViewport(0,0,w(),h());
        }
        glClear(GL_COLOR_BUFFER_BIT);

        xc = NV_Ith_S(state, 0);
        yc = NV_Ith_S(state, 1);
        theta = NV_Ith_S(state, 2);

        ox = 0.0;
        oy = 0.0;
        if (center == CENTER_OF_MASS) {
            ox = xc;
            oy = yc;
        }
        glPushMatrix();
            glColor3f(0.2, 0.8, 0.2);
            glBegin(GL_POLYGON);
            for (idx = 0; idx < 7; ++idx) {
                /* idx'th corner of the hexagon. */
                thetai = idx*M_PI/3 + theta;
                dx = params.L * cos(thetai);
                dy = params.L * sin(thetai);
                xx = (xc + dx - ox)/scale;
                yy = (yc + dy - oy)/scale;
                glVertex2f(xx, yy);
            }
            glEnd();

            /* Draw the disk in the center. */
            glColor3f(0.75, 0.75, 0.0);
            glBegin(GL_POLYGON);
            for (int i = 0; i < 13; i++) {
                double theta = 2*M_PI * i / 13;
                xx = (params.r0*cos(theta) - ox)/scale;
                yy = (params.r0*sin(theta) - oy)/scale;
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

        //std::cerr << "handle: e = " << e << "\n";
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
            else if (c == 'z') {
                scale *= 2.0;
            }
            else if (c == 'Z') {
                scale *= 0.5;
            }
        }
        return(Fl_Gl_Window::handle(e));
    }

    // Constructor
    Playback(int X, int Y, int W, int H, const char*L=0) : Fl_Gl_Window(X,Y,W,H,L)
    {
        int flag;
        double *p;

        frame = 0;
        frame_period = 1.0/24.0;

        center = ORIGIN;

        params.L = 1.25;
        params.g = 8.0;
        params.r0 = 0.5;

        state = N_VNew_Serial(6);
        NV_Ith_S(state, 0) = 9.0;
        NV_Ith_S(state, 1) = 9.0;
        NV_Ith_S(state, 2) = 0.0;

        NV_Ith_S(state, 3) =  0.45;
        NV_Ith_S(state, 4) = -0.45;
        NV_Ith_S(state, 5) =  0.0;


        tau = RCONST(0.0);
        dtau = RCONST(0.125);
        tau1 = RCONST(2500.0);

        cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
        if (cvode_mem == NULL) {
            std::cerr << "darn it\n";
        }

        flag = CVodeInit(cvode_mem, de_rigid_hex, tau, state);
        flag = CVodeSStolerances(cvode_mem, 1e-10, 1e-12);
        flag = CVodeSetUserData(cvode_mem, &params);
        flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
        flag = CVDense(cvode_mem, 6);
        flag = CVodeSetStopTime(cvode_mem, tau1);
        // flag = CVodeRootInit(cvode_mem, params.num_points, collision);

        crashed = 0;
        scale = 12.5;

        Fl::add_timeout(frame_period, Timer_CB, (void*)this);
        end();
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
