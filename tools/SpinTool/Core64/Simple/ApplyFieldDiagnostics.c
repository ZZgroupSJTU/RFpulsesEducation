#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: [Mx, My, Mz] = ApplyFieldDiagnostics(Bx, By, Bz, Min, tp)
//
// Simulates and returns the time evolution of a magnetization vector Min
// under the effect of a given field Bx(t), By(t), Bz(t) (Bi's are all
// vectors with N elements). If the field has N steps, then M will have 
// N+1 steps (MATLAB indicing used in diagram):
//
//
//         B(1)            B(2)                             B(N)
// M(1) ---------> M(2) ---------> M(3) -- ... ---> M(N) ---------> M(N+1)
// /|\                   dt=tp/N                                    /|\
//  |                                                                |
// t=0                                                              t=tp
// Input
//
// Input:
//      Field           Type            Units
//      Bx, By, Bz      1xN double      kHz
//      Min             3x1 double      a.u.
//      tp              1x1 double      ms
//
// Output:
//      Mx              1x(N+1)         
//      My              1x(N+1)         
//      Mz              1x(N+1)         
//
// IMPORTANT NOTE: RFamp must be real, BUT can also be negative (corresponding to phase = pi)
{
    int i,current_spin;                                    // counter variables
    const double pi = 3.1415926535897932384626433832795;   // Well, you know what this is
    double dt;                                             // Pulse time step in ms
    double tp;                                             // RF pulse duration
    int    NSteps;                                         // Number of pulse-steps (in time)
    double gm, cs;
    double *Min;
    double *Bx, *By, *Bz;
    double norm, angle, ct, st, nx, ny, nz, dotproduct;    // Internal variables for bloch equations
    double *Mxt, *Myt, *Mzt;
    
    
    // Set Gyromagnetic Ratio - we assume protons.
    // Units are kHz/Gauss/10
    gm = 4.257/10;
    
    // Retrieve spin data
    Bx  = mxGetPr(prhs[0]);
    By  = mxGetPr(prhs[1]);
    Bz  = mxGetPr(prhs[2]);
    Min = mxGetPr(prhs[3]);
    tp  = *mxGetPr(prhs[4]);
    
    // Number of (time-)steps in pulse
    NSteps = mxGetN(prhs[0]);
    
    // Time step, in ms
    dt = tp/NSteps;
    
    // Create output structure
    plhs[0] = mxCreateDoubleMatrix(1,NSteps+1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,NSteps+1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,NSteps+1,mxREAL);
    Mxt     = mxGetPr(plhs[0]);
    Myt     = mxGetPr(plhs[1]);
    Mzt     = mxGetPr(plhs[2]);
    Mxt[0]  = Min[0];
    Myt[0]  = Min[1];
    Mzt[0]  = Min[2];
    
    for (i=0; i<NSteps; i++)                                       // Loop over time (execute RF)
    {
        // calculate cosine and sine of rotation angles
        norm = sqrt(Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i]);      // Norm of B(t), in kHz
        if (norm==0)
            norm=1;                                                // Checks to see if we have division by 0, prolongs simulation by about 10%
        angle = -2*pi*norm*dt;                                     // Angle of rotation, in radians. Left hand rule.
        ct = cos(angle);
        st = sin(angle);
        nx = Bx[i]/norm;
        ny = By[i]/norm;
        nz = Bz[i]/norm;
        dotproduct = nx*Mxt[i] + ny*Myt[i] + nz*Mzt[i];
        // Use Rodriguez Formula for rotation:
        Mxt[i+1] = Mxt[i]*ct + (ny*Mzt[i] - nz*Myt[i])*st + dotproduct*(1-ct)*nx;
        Myt[i+1] = Myt[i]*ct + (nz*Mxt[i] - nx*Mzt[i])*st + dotproduct*(1-ct)*ny;
        Mzt[i+1] = Mzt[i]*ct + (nx*Myt[i] - ny*Mxt[i])*st + dotproduct*(1-ct)*nz;
    }

}