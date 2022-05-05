#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: [spinsOut, fid] = AcquirePulseRelax(spinsIn, pulse)
//
// Description: applies a given pulse to a spin ensemble, one spin at a time,
// while continuously acquiring the free induction decay (fid) from the
// sample. In particular, if the input pulse is null (having all entries = 0),
// this function can be used to simply acquire an FID. 
//
// This routine effectively solves the Bloch equations with relaxation, 
// under the simplifying assumption that the decay commutes with the rotation
// (valid for short time steps, or no RF).
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
//      1. A spin structure at the end of the pulse.
//      2. A complex FID, acquired at each point during the pulse.
//
//  [1] RFamp must be real, BUT can also be negative (corresponding to
//      a pi phase shift).
//  [2] This multiplies the FID acquired at the point.
{
    int i,current_spin;                                    // counter variables
    const double pi = 3.1415926535897932384626433832795;   // Well, you know what this is
    double vxR, vyR, vzR, vxI, vyI, vzI;                   // temporary vector components
    double dt;                                             // Pulse time step in ms
    double *RFamp, *RFphase, *Gx, *Gy, *Gz;                // pointers to RF data
    double tp;                                             // RF pulse duration
    int NSpins, NSteps;                                    // Number of spins and of pulse-steps
    double gm, cs, T1, T2, M0;
    double B0, B1, D; 
    double *r;
    double *Mout, *Min, *MoutIm, *MinIm;
    double *Bx, *By, *Bz, *fidRe, *fidIm;
    double norm, angle, ct, st, nx, ny, nz, dotproduct, R1, R2;    // Internal variables for bloch equations
    double RS;
    int    isComplex;
    
    // Set Gyromagnetic Ratio - we assume protons.
    // Units are kHz/Gauss/10
    gm = 4.257/10;
    
    // Retrieve Pulse data
    tp = *mxGetPr(mxGetField(prhs[1],0,"tp"));
    RFamp = mxGetPr(mxGetField(prhs[1],0,"RFamp"));
    RFphase = mxGetPr(mxGetField(prhs[1],0,"RFphase"));
    Gx = mxGetPr(mxGetField(prhs[1],0,"Gx"));
    Gy = mxGetPr(mxGetField(prhs[1],0,"Gy"));
    Gz = mxGetPr(mxGetField(prhs[1],0,"Gz"));
    
    NSteps = mxGetN(mxGetField(prhs[1],0,"RFamp"));
    // Time step, in ms
    dt = tp/NSteps;
    
    // Initialize arrays to contain RF x, y, z components (in kHz)
    Bx = (double*) malloc (sizeof(double)*NSteps);
    By = (double*) malloc (sizeof(double)*NSteps);
    Bz = (double*) malloc (sizeof(double)*NSteps);
    for (i=0; i<NSteps; i++)
    {
        Bx[i]=RFamp[i]*cos(RFphase[i]);
        By[i]=RFamp[i]*sin(RFphase[i]);
    }

    // Create output structure (identical to input structure)
    plhs[0] = mxDuplicateArray(prhs[0]);
    
    
    // Create output FID
    plhs[1] = mxCreateDoubleMatrix(1,NSteps,mxCOMPLEX);
    fidRe = mxGetPr(plhs[1]);
    fidIm = mxGetPi(plhs[1]);
    

    // Retrieve number of spins
    NSpins = mxGetNumberOfElements(prhs[0]);
    
    for (i=0; i<NSteps; i++)
    {
        fidRe[i] = 0;
        fidIm[i] = 0;
        for (current_spin=0; current_spin<NSpins; current_spin++)                   // Loop over spins
        {
            // Retrieve parameters for current spin & timestep
            Mout = mxGetPr(mxGetField(plhs[0], current_spin, "M"));
            if (mxIsComplex(mxGetField(prhs[0], current_spin, "M"))==1) {
                // It is possible for the input magnetization to be a complex
                // number, perhaps as a result of a phase cycling scheme
                // in which different cycles are added with complex coefficients.
                // However, to save on time, we only propagate the complex part
                // if it's actually there.
                MoutIm = mxGetPi(mxGetField(plhs[0], current_spin, "M")); 
                vxI = MoutIm[0]; vyI = MoutIm[1]; vzI = MoutIm[2];
                isComplex = 1;
            } else {
                isComplex = 0;
            }
            RS = mxGetScalar(mxGetField(prhs[0], current_spin, "RS")); // Scales FID
            B0 = mxGetScalar(mxGetField(prhs[0], current_spin, "B0")); // kHz
            B1 = mxGetScalar(mxGetField(prhs[0], current_spin, "B1")); // Scales RF
            // Sum Mxy
            fidRe[i] = fidRe[i] + Mout[0]*RS;
            fidIm[i] = fidIm[i] + Mout[1]*RS;
            if (isComplex==1) {
                // If the magnetization object is complex,
                //   Mx = MxR + i*MxI
                //   My = MyR + i*MyI
                // Then
                //   Mxy = Mx + i*My
                //       = (MxR - MyI) + i*(MxI + MyR)
                fidRe[i] = fidRe[i] - MoutIm[1]*RS;
                fidIm[i] = fidIm[i] + MoutIm[0]*RS;
            }
            vxR = Mout[0]; vyR = Mout[1]; vzR = Mout[2];
            r  = mxGetPr(mxGetField(prhs[0], current_spin, "r"));      // Get pointer to vector of spin location
            cs = mxGetScalar(mxGetField(prhs[0], current_spin, "cs")); // Get chemical shift of current spin (kHz)
            T1 = mxGetScalar(mxGetField(prhs[0], current_spin, "T1"));
            T2 = mxGetScalar(mxGetField(prhs[0], current_spin, "T2"));
            M0 = mxGetScalar(mxGetField(prhs[0], current_spin, "M0"));
            R2 = exp(-dt/T2);
            R1 = exp(-dt/T1);
            // Compute relevant parameters for current spin & timestep
            Bz[i] = cs + Gx[i]*r[0] + Gy[i]*r[1] + Gz[i]*r[2] + B0;    // in kHz

            norm = sqrt(B1*B1*Bx[i]*Bx[i] + B1*B1*By[i]*By[i] + Bz[i]*Bz[i]);
            angle = -2*pi*norm*dt;                                     // Angle of rotation, in radians. Note the minus sign - for the left-hand rule
            if (norm==0)
                norm=1;                                                // Checks to see if we have division by 0, elongates simulation by about 10%
            ct = cos(angle);
            st = sin(angle);
            nx = B1*Bx[i]/norm;
            ny = B1*By[i]/norm;
            nz = Bz[i]/norm;
            dotproduct = nx*vxR + ny*vyR + nz*vzR;
            // Use Rodriguez Formula for rotation:
            Mout[0] = vxR*ct + (ny*vzR - nz*vyR)*st + dotproduct*(1-ct)*nx;
            Mout[1] = vyR*ct + (nz*vxR - nx*vzR)*st + dotproduct*(1-ct)*ny;
            Mout[2] = vzR*ct + (nx*vyR - ny*vxR)*st + dotproduct*(1-ct)*nz;
            Mout[0] = Mout[0]*R2;
            Mout[1] = Mout[1]*R2;
            Mout[2] = (Mout[2]-M0)*R1 + M0;
            if (isComplex==1) {
                dotproduct = nx*vxI + ny*vyI + nz*vzI;
                // Use Rodriguez Formula for rotation:
                MoutIm[0] = vxI*ct + (ny*vzI - nz*vyI)*st + dotproduct*(1-ct)*nx;
                MoutIm[1] = vyI*ct + (nz*vxI - nx*vzI)*st + dotproduct*(1-ct)*ny;
                MoutIm[2] = vzI*ct + (nx*vyI - ny*vxI)*st + dotproduct*(1-ct)*nz;
                MoutIm[0] = MoutIm[0]*R2;
                MoutIm[1] = MoutIm[1]*R2;
                MoutIm[2] = (MoutIm[2]-M0)*R1 + M0;
            }
        }
    }
        
    
    // Free memory
    free(Bx);
    free(By);
    free(Bz);
    
}