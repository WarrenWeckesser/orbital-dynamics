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

class Playback : public Fl_Gl_Window {

    int frame;

    params_t params;
    N_Vector state;

    void *cvode_mem;
    realtype tau;
    realtype dtau;
    realtype tau1;


    //
    // Draw the system.
    //
    void draw()
    {
        double xx, yy;

        if (!valid()) {
            valid(1);
            glLoadIdentity();
            glViewport(0,0,w(),h());
        }
        glClear(GL_COLOR_BUFFER_BIT);
        /*
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_LINE_SMOOTH);
        glEnable(GL_POLYGON_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
        */
        glPushMatrix();
            glColor3f(0.2, 0.8, 0.2);
            // Draw the triangle with vertices
            // (state[0], state[1]), (state[2], state[3]), (state[4], state[5]) 
            //glBegin(GL_LINE_STRIP);
            glBegin(GL_POLYGON);
            for (int i = 0; i < 3; i++) {
                xx = X(state, i % 3)/12.5;
                yy = Y(state, i % 3)/12.5;
                glVertex2f(xx, yy);
            }
            glEnd();
            glColor3f(0.75, 0.75, 0.0);
            glBegin(GL_POLYGON);
            for (int i = 0; i < 13; i++) {
                double theta = 2*M_PI * i / 13;
                xx = 0.08*cos(theta);
                yy = 0.08*sin(theta);
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

        pb->redraw();
        Fl::repeat_timeout(1.0/60.0, Timer_CB, userdata);
    }

public:

    // Constructor
    Playback(int X, int Y, int W, int H, const char*L=0) : Fl_Gl_Window(X,Y,W,H,L)
    {
        int flag;

        frame = 0;

        params.k = 1.0;
        params.L = 1.0;
        params.b = 0.1;
        params.g = 9.0;

        state = N_VNew_Serial(12);
        X(state,0) = 8.9;
        Y(state,0)= 9.2;
        X(state,1) = 8.9;
        Y(state,1) = 8.2;
        X(state,2) = 9.3;
        Y(state,2) = 8.7;
        U(state,0) = 0.45;
        V(state,0) = -0.45;
        U(state,1) = 0.2;
        //V(state,1) = -0.25;
        V(state,1) = -2.5;
        U(state,2) = 0.5;
        V(state,2) = -0.4;


        tau = RCONST(0.0);
        dtau = RCONST(0.5);
        tau1 = RCONST(2500.0);

        cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
        if (cvode_mem == NULL) {
            std::cerr << "darn it\n";
        }

        flag = CVodeInit(cvode_mem, de3, tau, state);
        flag = CVodeSStolerances(cvode_mem, 1e-10, 1e-12);
        flag = CVodeSetUserData(cvode_mem, &params);
        flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
        flag = CVDense(cvode_mem, 12);
        flag = CVodeSetStopTime(cvode_mem, tau1);

        Fl::add_timeout(1.0/60.0, Timer_CB, (void*)this);
        end();
    }
};


int main()
{
    Fl_Window win(750, 750);
    Playback playback(10, 10, win.w()-20, win.h()-20);
    win.resizable(&playback);
    win.show();
    return(Fl::run());
}
