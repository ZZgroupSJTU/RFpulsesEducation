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
//      Field             Type                 Units
//      spins(i).r        3x1 double           mm
//      spins(i).M        3x1 double           a.u.
//      spins(i).cs       1x1 double           kHz
//      spins(i).T1       1x1 double           ms
//      spins(i).T2       1x1 double           ms
//      spins(i).M0       1x1 double           a.u.
//      spins(i).B1       1xNcT double/complex scaling (=1 for no change)
//      spins(i).B0       1x1 double           kHz (=offset)
//      spins(i).RS       1xNcR double/complex scaling (=1 for no change)
//
//      pulse structure (i=1,...,NcT):
//      Field             Type                 Units
//      pulse(i).tp       1x1 double           mm
//      pulse(i).RFamp    1xNp double          kHz [1]
//      pulse(i).RFphase  1xNp double          radians
//      pulse(i).Gx       1xNp double          kHz/mm
//      pulse(i).Gy       1xNp double          kHz/mm
//      pulse(i).Gz       1xNp double          kHz/mm
//
// Np - number of time points of RF
// NcT - number of transmit coils
// NcR - number of receive coils 
//
// Output:
//
//      A spin structure at the end of the pulse.
//
//  [1] RFamp must be real, BUT can also be negative (corresponding to
//      a pi phase shift).
{
    const double pi = 3.1415926535897932384626433832795;
    int numSpins, numSteps, numCoils;                   
    int curTimeIdx,curSpinIdx, curCoilIdx;
    double vx, vy, vz;                        // temporary vector components
    double dwellTime;                         // ms
    double *RFamp, *RFphase, *Gx, *Gy, *Gz;   // pointers to RF data
    double pulseDuration; 
    double cs;
    double T1, T2, B0, D; 
    double *B1real, *B1imag, *zeros;
    double B1Magnitude, B1Phase;
    double M0;                              
    double *r;
    double *Mout, *Min;
    double Bx, By, Bz;
    double effFieldMag, angle, ct, st, nx, ny, nz, dotproduct, R1, R2;    // Internal variables for bloch equations
    
    // Retrieve gradient data from the first pulse (they should all be the
    // same!)
    Gx = mxGetPr(mxGetField(prhs[1],0,"Gx"));
    Gy = mxGetPr(mxGetField(prhs[1],0,"Gy"));
    Gz = mxGetPr(mxGetField(prhs[1],0,"Gz"));
    
    // Retrieve number of steps and pulse duration using the first pulse
    pulseDuration = *mxGetPr(mxGetField(prhs[1],0,"tp")); // ms
    numSteps = mxGetN(mxGetField(prhs[1],0,"RFamp"));
    dwellTime = pulseDuration/numSteps;  // ms
    
    // Retrieve number of spins
    numSpins = mxGetNumberOfElements(prhs[0]);
    
    // Create output structure (identical to input structure)
    plhs[0] = mxDuplicateArray(prhs[0]);

    // Retrieve number of coils
    numCoils = mxGetNumberOfElements(prhs[1]);
    zeros = (double*) malloc (sizeof(double)*numCoils);
    for (curCoilIdx=0; curCoilIdx<numCoils; curCoilIdx++)
        zeros[curCoilIdx] = 0;
    
    
    for (curSpinIdx=0; curSpinIdx<numSpins; curSpinIdx++)          // Loop over spins
    {
        Mout = mxGetPr(mxGetField(plhs[0], curSpinIdx, "M"));  // Get pointer to magnetization of spin object
        Min  = mxGetPr(mxGetField(prhs[0], curSpinIdx, "M"));  // Get initial magnetization of current spin
        r    = mxGetPr(mxGetField(prhs[0], curSpinIdx, "r"));      // Get pointer to vector of spin location
        cs   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "cs")); // Get chemical shift of current spin (kHz)
        T1   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "T1"));
        T2   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "T2"));
        M0   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "M0"));
        B0   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "B0")); // kHz
        // RS   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "RS")); // kHz <-- not needed for transmission
        vx   = Min[0]; vy = Min[1]; vz = Min[2];
        if( mxIsComplex(mxGetField(prhs[0], curSpinIdx, "B1")) ) {
            B1real = mxGetPr(mxGetField(prhs[0], curSpinIdx, "B1")); // Scales RF
            B1imag = mxGetPi(mxGetField(prhs[0], curSpinIdx, "B1")); 
        } else {
            B1real = mxGetPr(mxGetField(prhs[0], curSpinIdx, "B1")); // Scales RF
            B1imag = zeros; 
        }
        
        // Calculate the field, as a function of time, felt by the spin, in kHz
        for (curTimeIdx=0; curTimeIdx<numSteps; curTimeIdx++)
        {
            // Calculate the instantaneous B field (in kHz) experienced by 
            // the current spin (as labeled by curSpinIdx)
            Bx = 0; By = 0;
            Bz = cs + B0 + Gx[curTimeIdx]*r[0] + Gy[curTimeIdx]*r[1] + Gz[curTimeIdx]*r[2];
            for (curCoilIdx=0; curCoilIdx<numCoils; curCoilIdx++)
            {
                RFamp   = mxGetPr(mxGetField(prhs[1],curCoilIdx,"RFamp"));
                RFphase = mxGetPr(mxGetField(prhs[1],curCoilIdx,"RFphase"));
                B1Magnitude = sqrt(B1real[curCoilIdx]*B1real[curCoilIdx] + B1imag[curCoilIdx]*B1imag[curCoilIdx]);
                B1Phase = atan2(B1imag[curCoilIdx], B1real[curCoilIdx]);
                Bx = Bx + B1Magnitude*RFamp[curTimeIdx]*cos(RFphase[curTimeIdx] + B1Phase);
                By = By + B1Magnitude*RFamp[curTimeIdx]*sin(RFphase[curTimeIdx] + B1Phase);
            }
            
            // calculate cosine and sine of rotation angles
            effFieldMag = sqrt(Bx*Bx + By*By + Bz*Bz);
            angle = -2*pi*effFieldMag*dwellTime;                                     // Angle of rotation, in radians. Note the minus sign - for the left-hand rule
            if (effFieldMag==0)
                effFieldMag=1;                                                // Checks to see if we have division by 0, elongates simulation by about 10%
            ct = cos(angle);
            st = sin(angle);
            // Compute the components of the instantaneous (and normalized) 
            // rotation axis
            nx = Bx/effFieldMag;
            ny = By/effFieldMag;
            nz = Bz/effFieldMag;
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
 
    free(zeros);
}