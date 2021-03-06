#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: spinsOut = ApplyPulse(spinsIn, pulse)
//
// Description: applies a given pulse to a spin ensemble, one spin at a time.
// This routine effectively solves the Bloch equations with relaxation, 
// under the simplifying assumption that the decay commutes with the rotation.
//
// Input:
//
//      spins structure:
//      Field           Type            Units
//      spins(i).r      3x1 double      mm
//      spins(i).M      3x1 double      a.u.
//      spins(i).cs     1x1 double      kHz
//      spins(i).T1     1x1 double      ms
//      spins(i).T2     1x1 double      ms
//      spins(i).M0     1x1 double      a.u.
//      spins(i).B1     1x1 double      scaling (=1 for no change)
//      spins(i).B0     1x1 double      kHz (=offset)
//
//      pulse structure:
//      Field           Type            Units
//      pulse.tp        1x1 double      mm
//      pulse.RFamp     1xNp double     kHz [1]
//      pulse.RFphase   1xNp double     radians
//      pulse.Gx        1xNp double     kHz/mm
//      pulse.Gy        1xNp double     kHz/mm
//      pulse.Gz        1xNp double     kHz/mm
//
// Output:
//
//      A spin structure at the end of the pulse.
//
//  [1] RFamp must be real, BUT can also be negative (corresponding to
//      a pi phase shift).
{
    const double pi = 3.1415926535897932384626433832795;
    int numSpins, numSteps;                   
    int curStepIdx,curSpinIdx;
    double vx, vy, vz;                        // temporary vector components
    double dwellTime;                         // ms
    double *RFamp, *RFphase, *Gx, *Gy, *Gz;   // pointers to RF data
    double pulseDuration; 
    double cs;
    double T1, T2, B0, B1, D; 
    double M0;                              
    double *r;
    double *Mout, *Min;
    double *Bx, *By, *Bz;
    double effFieldMag, angle, ct, st, nx, ny, nz, dotproduct, R1, R2;    // Internal variables for bloch equations
    
    // Retrieve Pulse data
    pulseDuration = *mxGetPr(mxGetFieldByNumber(prhs[1],0,0));
    RFamp = mxGetPr(mxGetFieldByNumber(prhs[1],0,1));
    RFphase = mxGetPr(mxGetFieldByNumber(prhs[1],0,2));
    Gx = mxGetPr(mxGetFieldByNumber(prhs[1],0,3));
    Gy = mxGetPr(mxGetFieldByNumber(prhs[1],0,4));
    Gz = mxGetPr(mxGetFieldByNumber(prhs[1],0,5));
    
    numSteps = mxGetN(mxGetFieldByNumber(prhs[1],0,1));
    dwellTime = pulseDuration/numSteps;  // ms
    
    // Initialize arrays to contain RF x, y, z components (in kHz)
    Bx = (double*) malloc (sizeof(double)*numSteps);
    By = (double*) malloc (sizeof(double)*numSteps);
    Bz = (double*) malloc (sizeof(double)*numSteps);
    for (curStepIdx=0; curStepIdx<numSteps; curStepIdx++)
    {
        Bx[curStepIdx]=RFamp[curStepIdx]*cos(RFphase[curStepIdx]);
        By[curStepIdx]=RFamp[curStepIdx]*sin(RFphase[curStepIdx]);
    }
    
    // Retrieve number of spins
    numSpins = mxGetNumberOfElements(prhs[0]);
    
    // Create output structure (identical to input structure)
    plhs[0] = mxDuplicateArray(prhs[0]);
    
    for (curSpinIdx=0; curSpinIdx<numSpins; curSpinIdx++)          // Loop over spins
    {
        Mout = mxGetPr(mxGetFieldByNumber(plhs[0], curSpinIdx, 1));  // Get pointer to magnetization of spin object
        Min = mxGetPr(mxGetFieldByNumber(prhs[0], curSpinIdx, 1));  // Get initial magnetization of current spin
        vx = Min[0]; vy = Min[1]; vz = Min[2];
        r = mxGetPr(mxGetFieldByNumber(prhs[0], curSpinIdx, 0));      // Get pointer to vector of spin location
        cs = mxGetScalar(mxGetFieldByNumber(prhs[0], curSpinIdx, 2)); // Get chemical shift of current spin (kHz)
        T1 = mxGetScalar(mxGetFieldByNumber(prhs[0], curSpinIdx, 3));
        T2 = mxGetScalar(mxGetFieldByNumber(prhs[0], curSpinIdx, 4));
        M0 = mxGetScalar(mxGetFieldByNumber(prhs[0], curSpinIdx, 5));
        B1 = mxGetScalar(mxGetFieldByNumber(prhs[0], curSpinIdx, 6)); // Scales RF
        B0 = mxGetScalar(mxGetFieldByNumber(prhs[0], curSpinIdx, 7)); // kHz
        for (curStepIdx=0; curStepIdx<numSteps; curStepIdx++)                                        // Calculate Bz in kHz as a function of time 
            Bz[curStepIdx] = cs + Gx[curStepIdx]*r[0] 
                                + Gy[curStepIdx]*r[1] 
                                + Gz[curStepIdx]*r[2]
                                + B0;
        for (curStepIdx=0; curStepIdx<numSteps; curStepIdx++)  // Execute RF
        {
            // calculate cosine and sine of rotation angles
            effFieldMag = sqrt(B1*B1*Bx[curStepIdx]*Bx[curStepIdx] + B1*B1*By[curStepIdx]*By[curStepIdx] + Bz[curStepIdx]*Bz[curStepIdx]);
            angle = -2*pi*effFieldMag*dwellTime;                                     // Angle of rotation, in radians. Note the minus sign - for the left-hand rule
            // angle = 2*pi*effFieldMag*dwellTime;                                     // Angle of rotation, in radians. Note the lack of minus sign - for the right-hand rule
            if (effFieldMag==0)
                effFieldMag=1;                                                // Checks to see if we have division by 0, elongates simulation by about 10%
            ct = cos(angle);
            st = sin(angle);
            // Compute the components of the instantaneous (and normalized) 
            // rotation axis
            nx = B1*Bx[curStepIdx]/effFieldMag;
            ny = B1*By[curStepIdx]/effFieldMag;
            nz = Bz[curStepIdx]/effFieldMag;
            dotproduct = nx*vx + ny*vy + nz*vz;
            // Use Rodriguez Formula for rotation:
            Mout[0] = vx*ct + (ny*vz - nz*vy)*st + dotproduct*(1-ct)*nx;
            Mout[1] = vy*ct + (nz*vx - nx*vz)*st + dotproduct*(1-ct)*ny;
            Mout[2] = vz*ct + (nx*vy - ny*vx)*st + dotproduct*(1-ct)*nz;
            R2 = exp(-dwellTime/T2);
            R1 = exp(-dwellTime/T1);
            Mout[0] = Mout[0]*R2;
            Mout[1] = Mout[1]*R2;
            Mout[2] = (Mout[2]-M0)*R1 + M0;
            vx = Mout[0];
            vy = Mout[1];
            vz = Mout[2];
        }
    }
    
    // Free memory
    free(Bx);
    free(By);
    free(Bz);
}