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
//      spins(i).RS     1x1 double      Receiver sensitivity, scaling (=1 for no change) [2]
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
//  [2] Only used when acquiring (not in this function)
{
    const double pi = 3.1415926535897932384626433832795;
    int numSpins, numSteps;                   
    int curStepIdx,curSpinIdx;
    double vxR, vyR, vzR, vxI, vyI, vzI;      // temporary vector components
    double dwellTime;                         // ms
    double *RFamp, *RFphase, *Gx, *Gy, *Gz;   // pointers to RF data
    double pulseDuration; 
    double cs;
    double T1, T2, B0, B1, D; 
    double M0;                              
    double *r;
    double *Mout, *Min, *MoutIm, *MinIm;
    double *Bx, *By, *Bz;
    double effFieldMag, angle, ct, st, nx, ny, nz, dotproduct, R1, R2;    // Internal variables for bloch equations
    double RS, B1amp, B1phase, B1real, B1imag, curBx, curBy;
    int    isComplex;                                                     // =1 if magnetization is complex. There are no boolean types in most C standards. Let's use int.
    int    is1DSim;                                                       // =1 if gradient is const, and all spins share the same T1, T2, B1, B0, RS
    int    isConstGrad;
    
    // Retrieve Pulse data
    pulseDuration = *mxGetPr(mxGetField(prhs[1],0,"tp"));
    RFamp = mxGetPr(mxGetField(prhs[1],0,"RFamp"));
    RFphase = mxGetPr(mxGetField(prhs[1],0,"RFphase"));
    Gx = mxGetPr(mxGetField(prhs[1],0,"Gx"));
    Gy = mxGetPr(mxGetField(prhs[1],0,"Gy"));
    Gz = mxGetPr(mxGetField(prhs[1],0,"Gz"));
    
    numSteps = mxGetN(mxGetField(prhs[1],0,"RFamp"));
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

    // Check for specialized case: if the gradient is constant, and all
    // the spins have the same T1, T2, B1, B0, RS, then a 1D simulation
    // can be carried out 
    /*
    is1DSim = 1;
    isConstGrad = 1;
    // First, check gradients are constant in time
    for (curStepIdx=0; curStepIdx<(numSteps-1); curStepIdx++)  {
        if (Gx[curStepIdx]!=Gx[curStepIdx+1]) { isConstGrad = 0; break; };
        if (Gy[curStepIdx]!=Gy[curStepIdx+1]) { isConstGrad = 0; break; };
        if (Gz[curStepIdx]!=Gz[curStepIdx+1]) { isConstGrad = 0; break; };
    }
    if (isConstGrad==1) {
        for (curSpinIdx=0; curSpinIdx<(numSpins-1); curSpinIdx++) {
            r    = mxGetPr(mxGetField(prhs[0], curSpinIdx, "r"));      // Get pointer to vector of spin location
            cs   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "cs")); // Get chemical shift of current spin (kHz)
            T1   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "T1"));
            T2   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "T2"));
            M0   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "M0"));
            if (mxIsComplex(mxGetField(prhs[0], curSpinIdx, "B1"))!=1)
            {
                B1amp = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "B1")); // Scales RF
                B1phase = 0.0; // Adds phase to RF
            } else {
                B1real = *mxGetPr(mxGetField(prhs[0], curSpinIdx, "B1")); 
                B1imag = *mxGetPi(mxGetField(prhs[0], curSpinIdx, "B1")); 
                B1amp = sqrt(B1real*B1real + B1imag*B1imag); // Scales RF
                B1phase = atan2(B1imag, B1real); // Adds phase to RF
            }
            B0 = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "B0")); // kHz
            RS = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "RS")); // kHz
        }
    }
    */
    
    for (curSpinIdx=0; curSpinIdx<numSpins; curSpinIdx++)          // Loop over spins
    {
        Mout = mxGetPr(mxGetField(plhs[0], curSpinIdx, "M"));      // Get pointer to magnetization of spin object
        // Min  = mxGetPr(mxGetField(prhs[0], curSpinIdx, "M"));      // Get initial magnetization of current spin
        vxR  = Mout[0]; vyR = Mout[1]; vzR = Mout[2];
        if (mxIsComplex(mxGetField(prhs[0], curSpinIdx, "M"))==1) {
            // It is possible for the input magnetization to be a complex
            // number, perhaps as a result of a phase cycling scheme
            // in which different cycles are added with complex coefficients.
            // However, to save on time, we only propagate the complex part
            // if it's actually there.
            // MinIm  = mxGetPi(mxGetField(prhs[0], curSpinIdx, "M")); 
            MoutIm = mxGetPi(mxGetField(plhs[0], curSpinIdx, "M")); 
            vxI = MoutIm[0]; vyI = MoutIm[1]; vzI = MoutIm[2];
            isComplex = 1;
        } else {
            isComplex = 0;
        }
        r    = mxGetPr(mxGetField(prhs[0], curSpinIdx, "r"));      // Get pointer to vector of spin location
        cs   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "cs")); // Get chemical shift of current spin (kHz)
        T1   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "T1"));
        T2   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "T2"));
        M0   = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "M0"));
        if (mxIsComplex(mxGetField(prhs[0], curSpinIdx, "B1"))!=1)
        {
            B1amp = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "B1")); // Scales RF
            B1phase = 0.0; // Adds phase to RF
        } else {
            B1real = *mxGetPr(mxGetField(prhs[0], curSpinIdx, "B1")); 
            B1imag = *mxGetPi(mxGetField(prhs[0], curSpinIdx, "B1")); 
            B1amp = sqrt(B1real*B1real + B1imag*B1imag); // Scales RF
            B1phase = atan2(B1imag, B1real); // Adds phase to RF
        }
        B0 = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "B0")); // kHz
        RS = mxGetScalar(mxGetField(prhs[0], curSpinIdx, "RS")); // kHz
        
        // Calculate Bz in kHz as a function of time 
        for (curStepIdx=0; curStepIdx<numSteps; curStepIdx++)                                        
            Bz[curStepIdx] = cs + Gx[curStepIdx]*r[0] 
                                + Gy[curStepIdx]*r[1] 
                                + Gz[curStepIdx]*r[2]
                                + B0;
        
        R2 = exp(-dwellTime/T2);
        R1 = exp(-dwellTime/T1);
        
        
        for (curStepIdx=0; curStepIdx<numSteps; curStepIdx++)  // Execute RF: Loop over time
        {
            // calculate cosine and sine of rotation angles
            if (B1phase!=0) {  // If B1 inho. has a phase, rotate transverse RF, RIGHT HAND RULE!
                curBx = Bx[curStepIdx]*cos(B1phase) - By[curStepIdx]*sin(B1phase);
                curBy = Bx[curStepIdx]*sin(B1phase) + By[curStepIdx]*sin(B1phase);
                curBx = curBx*B1amp;
                curBy = curBy*B1amp;
            } else {
                curBx = B1amp*Bx[curStepIdx];
                curBy = B1amp*By[curStepIdx];
            }
            effFieldMag = sqrt(curBx*curBx + curBy*curBy + Bz[curStepIdx]*Bz[curStepIdx]);
            angle = -2*pi*effFieldMag*dwellTime;                                 // Angle of rotation, in radians. Note the minus sign - for the left-hand rule
            if (effFieldMag==0) effFieldMag=1;                                   // Checks to see if we have division by 0, elongates simulation by about 10%
            ct = cos(angle);
            st = sin(angle);
            // Compute the components of the instantaneous (and normalized) 
            // rotation axis
            nx = curBx/effFieldMag;
            ny = curBy/effFieldMag;
            nz = Bz[curStepIdx]/effFieldMag;
            dotproduct = nx*vxR + ny*vyR + nz*vzR;
            // Use Rodriguez Formula for rotation:
            Mout[0] = vxR*ct + (ny*vzR - nz*vyR)*st + dotproduct*(1-ct)*nx;
            Mout[1] = vyR*ct + (nz*vxR - nx*vzR)*st + dotproduct*(1-ct)*ny;
            Mout[2] = vzR*ct + (nx*vyR - ny*vxR)*st + dotproduct*(1-ct)*nz;
            Mout[0] = Mout[0]*R2;
            Mout[1] = Mout[1]*R2;
            Mout[2] = (Mout[2]-M0)*R1 + M0;
            // Store results in temporary variable in preparation for next loop
            vxR = Mout[0];
            vyR = Mout[1];
            vzR = Mout[2];
            // If there is a complex component, do the exact same thing for it
            if (isComplex==1) {
                dotproduct = nx*vxI + ny*vyI + nz*vzI;
                // Use Rodriguez Formula for rotation:
                MoutIm[0] = vxI*ct + (ny*vzI - nz*vyI)*st + dotproduct*(1-ct)*nx;
                MoutIm[1] = vyI*ct + (nz*vxI - nx*vzI)*st + dotproduct*(1-ct)*ny;
                MoutIm[2] = vzI*ct + (nx*vyI - ny*vxI)*st + dotproduct*(1-ct)*nz;
                MoutIm[0] = MoutIm[0]*R2;
                MoutIm[1] = MoutIm[1]*R2;
                MoutIm[2] = (MoutIm[2]-M0)*R1 + M0;
                vxI = MoutIm[0];
                vyI = MoutIm[1];
                vzI = MoutIm[2];
            }
        }
    }
    
    // Free memory
    free(Bx);
    free(By);
    free(Bz);
}