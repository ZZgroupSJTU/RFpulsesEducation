#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: [Mx, My, Mz] = ApplyPulseDiagnosticsInverse(cs, r, Min, pulse)
// Evolves the Bloch equations BACKWARD in time, given a final state Min
// (The pulse will be played backwards) to yield an initial state. 
// Note that, if the pulse has N steps, then M will have N+1 steps
// (MATLAB indicing used in diagram):
//
//      pulse(1)         pulse(2)                        pulse(N)
// M(1) <--------- M(2) <--------- M(3) <- ... ---- M(N) <--------- M(N+1)
// /|\                   dt=tp/N                                    /|\
//  |                                                                |
// t=0                                                              t=tp
//                                                                  Input
//
//
//
// Input:
//      Field           Type            Units
//      cs              1x1 double      kHz
//      r               3x1 double      mm
//      Min             3x1 double      a.u.
// pulse structure:
//      pulse.tp        1x1 double      mm
//      pulse.RFamp     1xNp double     kHz
//      pulse.RFphase   1xNp double     radians
//      pulse.Gx        1xNp double     kHz/mm
//      pulse.Gy        1xNp double     kHz/mm
//      pulse.Gz        1xNp double     kHz/mm
//
// Output:
//      spinsout        spin structure, same as spins
//      Mt              3xNp            Matrix containing Mx(t), My(t), Mz(t)
//
// IMPORTANT NOTE: RFamp must be real, BUT can also be negative (corresponding to phase = pi)
{
    int i,current_spin;                                    // counter variables
    const double pi = 3.1415926535897932384626433832795;   // Well, you know what this is
    double vx, vy, vz;                                     // temporary vector components
    double dt;                                             // Pulse time step in ms
    double *RFamp, *RFphase, *Gx, *Gy, *Gz;                // pointers to RF data
    double tp;                                             // RF pulse duration
    int    NSteps;                                         // Number of pulse-steps (in time)
    double gm, cs;
    double *r;
    double *Min;
    double *Bx, *By, *Bz;
    double norm, angle, ct, st, nx, ny, nz, dotproduct;    // Internal variables for bloch equations
    double *Mxt, *Myt, *Mzt;
    
    
    // Set Gyromagnetic Ratio - we assume protons.
    // Units are kHz/Gauss/10
    gm = 4.257/10;
    
    // Retrieve spin data
    cs  = *mxGetPr(prhs[0]);
    r   = mxGetPr(prhs[1]);
    Min = mxGetPr(prhs[2]);
    
    // Retrieve Pulse data
    tp = *mxGetPr(mxGetFieldByNumber(prhs[3],0,0));
    RFamp = mxGetPr(mxGetFieldByNumber(prhs[3],0,1));
    RFphase = mxGetPr(mxGetFieldByNumber(prhs[3],0,2));
    Gx = mxGetPr(mxGetFieldByNumber(prhs[3],0,3));
    Gy = mxGetPr(mxGetFieldByNumber(prhs[3],0,4));
    Gz = mxGetPr(mxGetFieldByNumber(prhs[3],0,5));
    
    // Number of (time-)steps in pulse
    NSteps = mxGetN(mxGetFieldByNumber(prhs[3],0,1));
    
    // Time step, in ms
    dt = tp/NSteps;
    
    // Create output structure (identical to input structure)
    plhs[0] = mxCreateDoubleMatrix(1,NSteps+1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,NSteps+1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,NSteps+1,mxREAL);
    Mxt     = mxGetPr(plhs[0]);
    Myt     = mxGetPr(plhs[1]);
    Mzt     = mxGetPr(plhs[2]);
    Mxt[NSteps]  = Min[0];
    Myt[NSteps]  = Min[1];
    Mzt[NSteps]  = Min[2];
    
    // Initialize arrays to contain RF x, y, z components (in kHz)
    Bx = (double*) malloc (sizeof(double)*NSteps);
    By = (double*) malloc (sizeof(double)*NSteps);
    Bz = (double*) malloc (sizeof(double)*NSteps);
    for (i=0; i<NSteps; i++)
    {
        Bx[i]=RFamp[i]*cos(RFphase[i]);
        By[i]=RFamp[i]*sin(RFphase[i]);
        Bz[i]=cs + Gx[i]*r[0] + Gy[i]*r[1] + Gz[i]*r[2];
    }

    for (i=(NSteps-1); i>=0; i--)                                   // Loop BACK over time (execute RF)
    {
        // calculate cosine and sine of rotation angles
        norm = sqrt(Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i]);      // Norm of B(t), in kHz
        if (norm==0)
            norm=1;                                                // Checks to see if we have division by 0, prolongs simulation by about 10%
        angle = 2*pi*norm*dt;                                      // Angle of rotation, in radians. Inverted sign compare to the forward routine ApplyPulseDiagnostics!
        ct = cos(angle);
        st = sin(angle);
        nx = Bx[i]/norm;
        ny = By[i]/norm;
        nz = Bz[i]/norm;
        dotproduct = nx*Mxt[i+1] + ny*Myt[i+1] + nz*Mzt[i+1];
        // Use Rodriguez Formula for rotation:
        Mxt[i] = Mxt[i+1]*ct + (ny*Mzt[i+1] - nz*Myt[i+1])*st + dotproduct*(1-ct)*nx;
        Myt[i] = Myt[i+1]*ct + (nz*Mxt[i+1] - nx*Mzt[i+1])*st + dotproduct*(1-ct)*ny;
        Mzt[i] = Mzt[i+1]*ct + (nx*Myt[i+1] - ny*Mxt[i+1])*st + dotproduct*(1-ct)*nz;
    }
    
    // Free memory
    free(Bx);
    free(By);
    free(Bz);
}