#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: [spinsOut, fid] = ApplyPulse3D(spinsIn, pulse)
//
// Description: applies a given pulse to a 3D spin structure. 
// This routine effectively solves the Bloch equations with relaxation, 
// under the simplifying assumption that the decay commutes with the rotation.
//
// Input (some are shared by all spins):
//
//      spins structure:
//      Field           Type            Description
//      spins.B1      1x1 double        B1 Scaling (shared).
//      spins.cs     1x1 double         Chemical shift, kHz (shared).
//      spins.T1     1x1 double         T1 relaxation, ms (shared).
//      spins.T2     1x1 double         T2 relaxation, ms (shared).
//      spins.M0     1x1 double         Eq. magnetization, a.u. (shared).
//      spins.M      Nx*Ny*Nz*3 double  Magnetization vector at each pt.
//      spins.xVec   1xNx double        x-positions of spins in array, mm.
//      spins.yVec   1xNy double        y-positions of spins in array, mm.
//      spins.zVec   1xNz double        z-positions of spins in array, mm.
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
//      An fid, obtained by summing all x & y magnetization at each pulse step.
//
//  [1] RFamp must be real, BUT can also be negative (corresponding to
//      a pi phase shift).
//
// Theory:
// The simulation calculates a set of propagators (3x3 matrices), valid
// for a given offset, T1, T2, B1, etc ... value. Under the assumption of
// a short dwell time, the evolution is broken up into concatenated steps
// of rotation followed by relaxation. In terms of operators,
//
//          Rotation               Relax
//    M(i) ----------> R(i)*M(i) --------> E*R(i)*M(i) + B
//
// where R(i) is the 3x3 rotation matrix associated with the i-th step,
// and E is a constant relaxation matrix
//
//       [E2 0  0 ]        [0        ]
//   E = [0  E2 0 ]    B = [0        ]    E1 = exp(-dt/T1)   E2 = exp(-dt/T2)
//       [0  0  E1]        [(1-E1)*M0]
//
// The effect of several such concatenated steps is (denoting Q(i)=E*R(i)
//
//   M(i) ---> M(i+1) = E*R(i)*M(i) + B = Q(i)*M(i) + B
//        ---> M(i+2) = Q(i+1)*M(i+1) + B
//                    = Q(i+1)*(Q(i)*M(i) + B) + B
//                    = Q(i+1)*Q(i)*M(i) + (1+Q(i+1))*B
//        ---> M(i+3) = Q(i+2)*M(i+2) + B
//                    = Q(i+2)*(Q(i+1)*Q(i)*M(i) + (1+Q(i+1))*B) + B
//                    = Q(i+2)*Q(i+1)*E*R(i)*M(i) + Q(i+2)*(1+Q(i+1))*B + B
//                    = Q(i+2)*Q(i+1)*Q(i)*M(i) + (1 + Q(i+2) + Q(i+2)*Q(i+1))*B
//        ---> M(i+4) = Q(i+3)*Q(i+2)*Q(i+1)*Q(i)*M(i) + (1+Q(i+3)+Q(i+3)*Q(i+2)+Q(i+3)*Q(i+2)*Q(i+1))*B
//
// et cetera. Thus, this is a linear transformation of the form
//
//   M(i) ---> M(i+N) = F1*M(i) + F2
//   F1 = Q(N-1)*...*Q(i+3)*Q(i+2)*Q(i+1)*Q(i)
//   F2 = 1 + Q(N-1) + Q(N-1)*Q(N-2) + ... + Q(N-1)*Q(i+3)*Q(i+2)*Q(i+1)
//
// which can be applied to any initial magnetization M(i).
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
    int    numOffsets;
    double **offsetVec;
    double *xVec, *yVec, *zVec;
    int    idxOffset, idxTime;
    int    pulseAxis;   
    int    Nx, Ny, Nz;
    double *Mi, *Mo;
    double Q[3][3];
    double F1[3][3];
    double F2[3][3];
    double F1Temp[3][3];
    double E[3], ETemp[3];
    double T;
    long   iX, iY, iZ, prodX, prodXY, prodXYZ, idxX, idxY, idxZ, Mx, My, Mz, p;
    
    // Retrieve Pulse data
    pulseDuration = *mxGetPr(mxGetField(prhs[1],0,"tp"));
    RFamp = mxGetPr(mxGetField(prhs[1],0,"RFamp"));
    RFphase = mxGetPr(mxGetField(prhs[1],0,"RFphase"));
    Gx = mxGetPr(mxGetField(prhs[1],0,"Gx"));
    Gy = mxGetPr(mxGetField(prhs[1],0,"Gy"));
    Gz = mxGetPr(mxGetField(prhs[1],0,"Gz"));
    numSteps = mxGetN(mxGetField(prhs[1],0,"RFamp"));
    dwellTime = pulseDuration/numSteps;  // ms

    // Retrieve spin data
    T1 = *mxGetPr(mxGetField(prhs[0],0,"T1"));
    T2 = *mxGetPr(mxGetField(prhs[0],0,"T2"));
    B1 = *mxGetPr(mxGetField(prhs[0],0,"B1"));
    cs = *mxGetPr(mxGetField(prhs[0],0,"cs"));
    M0 = *mxGetPr(mxGetField(prhs[0],0,"M0"));
    xVec = mxGetPr(mxGetField(prhs[0],0,"xVec"));
    yVec = mxGetPr(mxGetField(prhs[0],0,"yVec"));
    zVec = mxGetPr(mxGetField(prhs[0],0,"zVec"));
    Nx = mxGetN(mxGetField(prhs[0],0,"xVec"));
    Ny = mxGetN(mxGetField(prhs[0],0,"yVec"));
    Nz = mxGetN(mxGetField(prhs[0],0,"zVec"));
    R2 = exp(-dwellTime/T2);
    R1 = exp(-dwellTime/T1);

   
    // Initialize arrays to contain RF x, y, z components (in kHz)
    Bx = (double*) malloc (sizeof(double)*numSteps);
    By = (double*) malloc (sizeof(double)*numSteps);
    for (curStepIdx=0; curStepIdx<numSteps; curStepIdx++)
    {
        Bx[curStepIdx]=B1*RFamp[curStepIdx]*cos(RFphase[curStepIdx]); // kHz
        By[curStepIdx]=B1*RFamp[curStepIdx]*sin(RFphase[curStepIdx]); // kHz
    }

    // Create output structure (identical to input structure)
    plhs[0] = mxDuplicateArray(prhs[0]);
    Mi = mxGetPr(mxGetField(prhs[0], 0, "M"));  
    Mo = mxGetPr(mxGetField(plhs[0], 0, "M"));  

    
    // Check gradients are constant in time
    for (curStepIdx=0; curStepIdx<(numSteps-1); curStepIdx++)  {
        if (Gx[curStepIdx]!=Gx[curStepIdx+1]) { isConstGrad = 0; break; };
        if (Gy[curStepIdx]!=Gy[curStepIdx+1]) { isConstGrad = 0; break; };
        if (Gz[curStepIdx]!=Gz[curStepIdx+1]) { isConstGrad = 0; break; };
    }
    
    pulseAxis = -1; 
    offsetVec = (double**) malloc (sizeof(double*)*numSteps);
    
    if ((Gx[0]!=0) & (Gy[0]==0) & (Gz[0]==0)) { 
        pulseAxis = 1; 
        numOffsets = Nx;
        for (idxTime=0; idxTime<numSteps; idxTime++) {
            offsetVec[idxTime] = (double*) malloc (sizeof(double)*Nx);
            for (idxOffset=0; idxOffset<numOffsets; idxOffset++) {
                offsetVec[idxTime][idxOffset] = cs + Gx[idxTime]*xVec[idxOffset];
            }
        }
    };
    
    if ((Gx[0]==0) & (Gy[0]!=0) & (Gz[0]==0)) { 
        pulseAxis = 2; 
        numOffsets = Ny;
        for (idxTime=0; idxTime<numSteps; idxTime++) {
            offsetVec[idxTime] = (double*) malloc (sizeof(double)*Ny);
            for (idxOffset=0; idxOffset<numOffsets; idxOffset++) {
                offsetVec[idxTime][idxOffset] = cs + Gy[idxTime]*yVec[idxOffset];
            }
        }
    };
    
    if ((Gx[0]==0) & (Gy[0]==0) & (Gz[0]!=0)) { 
        pulseAxis = 3; 
        numOffsets = Nz;
        for (idxTime=0; idxTime<numSteps; idxTime++) {
            offsetVec[idxTime] = (double*) malloc (sizeof(double)*Nz);
            for (idxOffset=0; idxOffset<numOffsets; idxOffset++) {
                offsetVec[idxTime][idxOffset] = cs + Gz[idxTime]*zVec[idxOffset];
            }
        }
    };
    
    if ((Gx[0]==0) & (Gy[0]==0) & (Gz[0]==0)) { 
        pulseAxis = 0; 
        numOffsets = 1;
        offsetVec[idxTime] = (double*) malloc(sizeof(double)*1);
        for (idxTime=0; idxTime<numSteps; idxTime++) {
            offsetVec[idxTime] = (double*) malloc (sizeof(double)*1);
            offsetVec[idxTime][0] = cs;
        }
    };
    
    T = (1-R1)*M0;
    if (pulseAxis>=0) {
        for (idxOffset=0; idxOffset<numOffsets; idxOffset++) {
            F1[0][0] = 1;  F1[0][1] = 0;  F1[0][2] = 0;
            F1[1][0] = 0;  F1[1][1] = 1;  F1[1][2] = 0;
            F1[2][0] = 0;  F1[2][1] = 0;  F1[2][2] = 1;
            F2[0][0] = 1;  F2[0][1] = 0;  F2[0][2] = 0;
            F2[1][0] = 0;  F2[1][1] = 1;  F2[1][2] = 0;
            F2[2][0] = 0;  F2[2][1] = 0;  F2[2][2] = 1;
            E[0] = 0;
            E[1] = 0;
            E[2] = T;
            // We calculate the "forward propagator" Mi-->F1*Mi + F2
            // This is best done *in reverse*, starting from the last step
            // and building up F1 and F2 "backward".
            for (idxTime=numSteps-1; idxTime>=0; idxTime--) {
                effFieldMag = sqrt(Bx[idxTime]*Bx[idxTime] + By[idxTime]*By[idxTime] + offsetVec[idxTime][idxOffset]*offsetVec[idxTime][idxOffset]);
                angle = -2*pi*effFieldMag*dwellTime;                                 // Angle of rotation, in radians. Note the minus sign - for the left-hand rule
                ct = cos(angle);
                st = sin(angle);
                // Compute the components of the instantaneous (and normalized) 
                // rotation axis
                if (effFieldMag==0) {
                    Q[0][0] = R2;    Q[0][1] = 0;    Q[0][2] = 0;
                    Q[1][0] = 0;     Q[1][1] = R2;   Q[1][2] = 0;
                    Q[2][0] = 0;     Q[2][1] = 0;    Q[2][2] = R1;
                } else {
                    nx = Bx[idxTime]/effFieldMag;
                    ny = By[idxTime]/effFieldMag;
                    nz = offsetVec[idxTime][idxOffset]/effFieldMag;
                    Q[0][0] = (ct+nx*nx*(1-ct))*R2;     Q[0][1] = (nx*ny*(1-ct)-nz*st)*R2;   Q[0][2] = (nx*nz*(1-ct)+ny*st)*R2;
                    Q[1][0] = (ny*nx*(1-ct)+nz*st)*R2;  Q[1][1] = (ct+ny*ny*(1-ct))*R2;      Q[1][2] = (ny*nz*(1-ct)-nx*st)*R2;
                    Q[2][0] = (nz*nx*(1-ct)-ny*st)*R1;  Q[2][1] = (nz*ny*(1-ct)+nx*st)*R1;   Q[2][2] = (ct+nz*nz*(1-ct))*R1;
                }
                
                F1Temp[0][0] = F1[0][0];   F1Temp[0][1] = F1[0][1];    F1Temp[0][2] = F1[0][2];    
                F1Temp[1][0] = F1[1][0];   F1Temp[1][1] = F1[1][1];    F1Temp[1][2] = F1[1][2];    
                F1Temp[2][0] = F1[2][0];   F1Temp[2][1] = F1[2][1];    F1Temp[2][2] = F1[2][2];
                // Multiply F1(i+1) = F1(i)*Q(i)
                F1[0][0] = F1Temp[0][0]*Q[0][0] + F1Temp[0][1]*Q[1][0] + F1Temp[0][2]*Q[2][0];
                F1[0][1] = F1Temp[0][0]*Q[0][1] + F1Temp[0][1]*Q[1][1] + F1Temp[0][2]*Q[2][1];
                F1[0][2] = F1Temp[0][0]*Q[0][2] + F1Temp[0][1]*Q[1][2] + F1Temp[0][2]*Q[2][2];
                F1[1][0] = F1Temp[1][0]*Q[0][0] + F1Temp[1][1]*Q[1][0] + F1Temp[1][2]*Q[2][0];
                F1[1][1] = F1Temp[1][0]*Q[0][1] + F1Temp[1][1]*Q[1][1] + F1Temp[1][2]*Q[2][1];
                F1[1][2] = F1Temp[1][0]*Q[0][2] + F1Temp[1][1]*Q[1][2] + F1Temp[1][2]*Q[2][2];
                F1[2][0] = F1Temp[2][0]*Q[0][0] + F1Temp[2][1]*Q[1][0] + F1Temp[2][2]*Q[2][0];
                F1[2][1] = F1Temp[2][0]*Q[0][1] + F1Temp[2][1]*Q[1][1] + F1Temp[2][2]*Q[2][1];
                F1[2][2] = F1Temp[2][0]*Q[0][2] + F1Temp[2][1]*Q[1][2] + F1Temp[2][2]*Q[2][2];
                // Add F2 => F2 + F1
                if (idxTime>0) {
                    E[0] = E[0] + F1[0][2]*T; 
                    E[1] = E[1] + F1[1][2]*T; 
                    E[2] = E[2] + F1[2][2]*T; 
                };
                /*
                mexPrintf("F1[0][0] = %.2f   F1[0][1] = %.2f   F1[0][2] = %.2f \n", F1[0][0], F1[0][1], F1[0][2]);
                mexPrintf("F1[1][0] = %.2f   F1[1][1] = %.2f   F1[1][2] = %.2f \n", F1[1][0], F1[1][1], F1[1][2]);
                mexPrintf("F1[2][0] = %.2f   F1[2][1] = %.2f   F1[2][2] = %.2f \n", F1[2][0], F1[2][1], F1[2][2]);
                 */
                
            }
            // Now that the propagator has been calculated, let us apply it to all spins having this particular offset
            // A note about MATLAB indexing: If a matrix M[i][j][k][p] is passed onto a mex file, with
            // dimension sizes Nx*Ny*Nz*Nt, then to access the (i,j,k,p) index we need to access
            //    (i,j,k,p) ---->  Linear index:  (i-1) + (j-1)*Nx + (k-1)*Nx*Ny + (p-1)*Nx*Ny*Nz
            prodX = Nx;
            prodXY = Nx*Ny;
            prodXYZ = Nx*Ny*Nz;
            switch (pulseAxis) {
                case 0:
                    for (idxX=0; idxX<Nx; idxX++) {
                        iX = idxX;
                        for (idxY=0; idxY<Ny; idxY++) {
                            iY = idxY*prodX;
                            for (idxZ=0; idxZ<Nz; idxZ++) {
                                iZ = idxZ*prodXY;
                                Mo[iX+iY+iZ+0]         = E[0] + F1[0][0]*Mi[iX+iY+iZ] + F1[0][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[0][2]*Mi[iX+iY+iZ+2*prodXYZ];
                                Mo[iX+iY+iZ+1*prodXYZ] = E[1] + F1[1][0]*Mi[iX+iY+iZ] + F1[1][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[1][2]*Mi[iX+iY+iZ+2*prodXYZ];
                                Mo[iX+iY+iZ+2*prodXYZ] = E[2] + F1[2][0]*Mi[iX+iY+iZ] + F1[2][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[2][2]*Mi[iX+iY+iZ+2*prodXYZ];
                            }
                        }
                    }
                    break;
                case 1:
                    // Pulse is along X
                    iX = idxOffset;
                    for (idxY=0; idxY<Ny; idxY++) {
                        iY = idxY*prodX;
                        for (idxZ=0; idxZ<Nz; idxZ++) {
                            iZ = idxZ*prodXY;
                            Mo[iX+iY+iZ+0]         = E[0] + F1[0][0]*Mi[iX+iY+iZ] + F1[0][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[0][2]*Mi[iX+iY+iZ+2*prodXYZ];
                            Mo[iX+iY+iZ+1*prodXYZ] = E[1] + F1[1][0]*Mi[iX+iY+iZ] + F1[1][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[1][2]*Mi[iX+iY+iZ+2*prodXYZ];
                            Mo[iX+iY+iZ+2*prodXYZ] = E[2] + F1[2][0]*Mi[iX+iY+iZ] + F1[2][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[2][2]*Mi[iX+iY+iZ+2*prodXYZ];
                        }
                    }
                    break;
                case 2:
                    // Pulse is along Y
                    iY = idxOffset*prodX;
                    for (idxX=0; idxX<Nx; idxX++) {
                        iX = idxX;
                        for (idxZ=0; idxZ<Nz; idxZ++) {
                            iZ = idxZ*prodXY;
                            Mo[iX+iY+iZ+0]         = E[0] + F1[0][0]*Mi[iX+iY+iZ] + F1[0][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[0][2]*Mi[iX+iY+iZ+2*prodXYZ];
                            Mo[iX+iY+iZ+1*prodXYZ] = E[1] + F1[1][0]*Mi[iX+iY+iZ] + F1[1][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[1][2]*Mi[iX+iY+iZ+2*prodXYZ];
                            Mo[iX+iY+iZ+2*prodXYZ] = E[2] + F1[2][0]*Mi[iX+iY+iZ] + F1[2][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[2][2]*Mi[iX+iY+iZ+2*prodXYZ];
                        }
                    }
                    break;
                case 3:
                    // Pulse is along Z
                    iZ = idxOffset*prodXY;
                    for (idxX=0; idxX<Nx; idxX++) {
                        iX = idxX;
                        for (idxY=0; idxY<Ny; idxY++) {
                            iY = idxY*prodX;
                            Mo[iX+iY+iZ+0]         = E[0] + F1[0][0]*Mi[iX+iY+iZ] + F1[0][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[0][2]*Mi[iX+iY+iZ+2*prodXYZ];
                            Mo[iX+iY+iZ+1*prodXYZ] = E[1] + F1[1][0]*Mi[iX+iY+iZ] + F1[1][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[1][2]*Mi[iX+iY+iZ+2*prodXYZ];
                            Mo[iX+iY+iZ+2*prodXYZ] = E[2] + F1[2][0]*Mi[iX+iY+iZ] + F1[2][1]*Mi[iX+iY+iZ+1*prodXYZ] + F1[2][2]*Mi[iX+iY+iZ+2*prodXYZ];
                        }
                    }
                    break;
            }
        }
        for (idxTime=0; idxTime<numSteps; idxTime++) {
            free(offsetVec[idxTime]);
        }
        free(offsetVec);
    };

    /*
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
    */
    
    
    // Free memory
    free(Bx);
    free(By);
}

int multInts(int x, int y) {
    return x*y;
}