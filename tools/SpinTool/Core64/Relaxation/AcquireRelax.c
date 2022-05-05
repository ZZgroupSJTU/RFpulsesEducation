#include "mex.h"
#include <math.h>
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: [spinsOut,fid] = AcquireRelax(spinsIn, pulse)
//
// Input:
//      Field           Type            Units
// spins structure:
//      spins(i).r      3x1 double      mm
//      spins(i).M      3x1 double      a.u.
//      spins(i).cs     1x1 double      kHz
//      spins(i).T1     1x1 double      ms
//      spins(i).T2     1x1 double      ms
//      spins(i).M0     1x1 double      a.u.
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
//      fid             2xNp            a.u.
//
// IMPORTANT NOTE: This function uses only the gradients within the pulse
// structure, and disregards the RF. Consequently, all spins are assumed to
// have Mz=constant and we work with 2D vectors.
{
    int i,current_spin;                                    // counter variables (i used for time)
    const double pi = 3.1415926535897932384626433832795;   // Well, you know what this is
    double dt;                                             // Pulse time step in ms
    double *Gx, *Gy, *Gz;                                  // pointers to gradient components
    double tp;                                             // RF pulse duration (=acquisition duration)
    int    NSpins, NSteps;                                 // Number of spins and of pulse-steps
    double cs;
    double *r;
    double *Mout, *Min;
    double *Phases, *Mxy;                                  // Phases & magnitudes of xy spins as func. of position
    double *FIDReal, *FIDImag;                             // Pointer to output FID. 
    double tempFID_r, tempFID_c;
    double dPhase;
    double T2, T1, M0;
    double expT1;
    
    // Retrieve Pulse data
    tp = *mxGetPr(mxGetFieldByNumber(prhs[1],0,0));
    Gx = mxGetPr(mxGetFieldByNumber(prhs[1],0,3));
    Gy = mxGetPr(mxGetFieldByNumber(prhs[1],0,4));
    Gz = mxGetPr(mxGetFieldByNumber(prhs[1],0,5));
    NSteps = mxGetN(mxGetFieldByNumber(prhs[1],0,1));
    // Time step, in ms
    dt = tp/NSteps;

    // Retrieve number of spins
    NSpins = mxGetNumberOfElements(prhs[0]);

    // Initialize array to contain phases of xy spins as func. of position (in radians)
    Phases = (double*) malloc (sizeof(double)*NSpins);
    Mxy    = (double*) malloc (sizeof(double)*NSpins);
    
    // Retrieve initial phases of spins (between 0 and 2*pi)
    for (current_spin=0; current_spin<NSpins; current_spin++)
    {
        Min = mxGetPr(mxGetFieldByNumber(prhs[0], current_spin, 1));  // Get pointer to magnetization of spin object
        Phases[current_spin] = atan2(Min[1],Min[0]);
        if (Phases[current_spin]<0) 
           Phases[current_spin] = Phases[current_spin];
        Mxy[current_spin] = sqrt(Min[0]*Min[0] + Min[1]*Min[1]);
    }

    // Create output spin structure (identical to input structure) & FID matrix (initialized to 0)
    plhs[0] = mxDuplicateArray(prhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,NSteps,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,NSteps,mxREAL);
    FIDReal = mxGetPr(plhs[1]);
    FIDImag = mxGetPr(plhs[2]);

    // Acquire
    for (i=0; i<NSteps; i++)                                                        // Loop over time
    {
        // sum Mxy
        tempFID_r = 0;
        tempFID_c = 0;
        for (current_spin=0; current_spin<NSpins; current_spin++)                   // Loop over spins
        {
            tempFID_r = tempFID_r + Mxy[current_spin]*cos(Phases[current_spin]);
            tempFID_c = tempFID_c + Mxy[current_spin]*sin(Phases[current_spin]);
        }
        FIDReal[i] = tempFID_r;
        FIDImag[i] = tempFID_c;
        // Propagate spins 
        for (current_spin=0; current_spin<NSpins; current_spin++)                   // Loop over spins
        {
            // Get parameters of current spin
            r  = mxGetPr(mxGetFieldByNumber(prhs[0], current_spin, 0));              // Get pointer to vector of spin location (mm)
            cs = mxGetScalar(mxGetFieldByNumber(prhs[0], current_spin, 2));          // Get chemical shift of current spin (kHz)
            T1 = mxGetScalar(mxGetFieldByNumber(prhs[0], current_spin, 3));
            T2 = mxGetScalar(mxGetFieldByNumber(prhs[0], current_spin, 4));
            M0 = mxGetScalar(mxGetFieldByNumber(prhs[0], current_spin, 5));
            // increment Mxy phase and decrease Mxy magnitude by e^(-t/T2)
            // Note dPhase has a minus sign to keep the left-hand convention of spins
            dPhase = -2*pi*dt*(cs + Gx[i]*r[0] + Gy[i]*r[1] + Gz[i]*r[2]); 
            Phases[current_spin] = Phases[current_spin] + dPhase;
            Mxy[current_spin] = Mxy[current_spin]*exp(-dt/T2);
            // Relax Mz towards M0 (Mz(t+dt) = e^(-t/T1)*Mz(t) + (1-e^(-t/T1))*M0)
            Mout = mxGetPr(mxGetFieldByNumber(plhs[0], current_spin, 1));
            expT1 = exp(-dt/T1);
            Mout[2] = Mout[2]*expT1 + (1-expT1)*M0;
        }
    }

    // Set output spins properly
    for (current_spin=0; current_spin<NSpins; current_spin++)
    {
        Mout    = mxGetPr(mxGetFieldByNumber(plhs[0], current_spin, 1));
        Mout[0] = Mxy[current_spin]*cos(Phases[current_spin]); 
        Mout[1] = Mxy[current_spin]*sin(Phases[current_spin]); 
    }
    
    // Free memory
    free(Phases);
    free(Mxy);
}