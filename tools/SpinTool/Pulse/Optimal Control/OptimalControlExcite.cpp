#include "mex.h"
#include <math.h>
#include "armadillo"

using namespace arma;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: RF = OptimalControl(RF0, tp, offsetVec, B1Vec, targetU, maxIter, epsilon, isOutput, maxB1)
//
//      Variable Name   Type            Units     Description
//
// Input:
//      RF0        1xN complex          kHz       Initial RF guess
//      tp         1x1 double           ms        Pulse duration (fixed)
//      offsetVec  1xM double           kHz       offsets to optimize for
//      B1Vec      1xP double           [scale]   B1 values to optimize for
//      targetRho  {PxM}x(2x2) complex  -         Target propagator, as func. of offset & B1
//      maxIter    1x1 int              -         Max. number of iterations
//      epsilon    1x1 double           -         Step size, should be >0 and small
//      isOutput   1x1 boolean (0, 1)   -         Output feedback (if set to 1)
//      maxB1      1x1 double           kHz       Max. RF amplitude
//
// Output:
//      RF              1xN            Matrix containing Mx(t), My(t), Mz(t)
//
// Notes:
//
// Forward propagation: X(j) = U(j)*X(j-1)
//
// t = 0        1           2                  N-1
//     U1       U2          U3                 UN
//   I ----> U1 ----> U2*U1 ----> U3*U2*U1 ... ----> UN*...*U1
//   X0      X1       X2=U2*X1    X3=U3*X2           XN=UN*X(N-1)
//
// Backward propagation: P(j) = U(j)*P(j-1)
// ! = dagger 
//               
//   I  <---- UN!    ...     <--- U3!*...*UN!  <--- U2!*...*UN!  <--- U1!*U2!*...*UN! 
//   PN       P(N-1)              P2=U2*P1          P1=U1*P0          P0=XN!*UF
//             =U(N-1)*P(N-2)                                                        
//
// 
// Connection with notation in Glaser, 2005:
// Lj = Lambda_j = Pj*C*Xj
// Rj = Rho_j    = 
{
    const double pi = 3.1415926535897932384626433832795;   // Well, you know what this is
    int b; // B1 index
    int f; // Offset index
    int t; // Time index
    int isDebug = 0;
    int idxIter = 0;
    int isConverged = 0;
    int isOutput = 1;
    int maxIter = 1000;
    int numB1s = 5;
    int numOffsets = 50;
    int N = 1; // Elements in RF
    int i;
    double nx, ny, nz; // Rotation axis, nx^2 + ny^2 + nz^2 = 1
    double phi; // Rotation angle
    double Bx, By, Bz;
    double gamma = 42.576; // kHz/mT = MHz/T
    double dt = 1; // Dwell time, in ms
    double *RFX0, *RFY0;
    double epsilon;
    double aRe, aIm, bRe, bIm;
    double sp, cp; // sin(phi/2), cos(phi/2)
    double a1, a2, a3, a4, a5, a6, q1;
    double tp;
    double *offsetVec, *B1Vec;
    double *UTRe;
    double *UTIm;
    double *RFX, *RFY;
    double Beff;
    double maxB1, magB1;
    double curCost, prevCost;
    double twoZeros[2];
    int subs0[2];
    mwSize ndims = 2;
    
    cx_mat U(2,2);
    cx_mat targetU(2,1);
    double c1, c2, c3;
    cx_mat fullX(2,2);
    cx_mat fullP(2,2);
    
    // Initialize auxiliary structures, of size (N+1)*2
//     double **XRe;
//     double **XIm;
//     double **PRe;
//     double **PIm;

    twoZeros[0] = 0;
    twoZeros[1] = 0;
    
    // Retrieve inputs
    N = mxGetN(prhs[0]);
    RFX0 = mxGetPr(prhs[0]);
    RFY0 = mxGetPi(prhs[0]);
    tp = *mxGetPr(prhs[1]);
    offsetVec = mxGetPr(prhs[2]);
    B1Vec = mxGetPr(prhs[3]);
    maxIter = (int)(*mxGetPr(prhs[5])); 
    epsilon = *mxGetPr(prhs[6]); 
    isOutput = *mxGetPr(prhs[7]); 
    maxB1 = *mxGetPr(prhs[8]);

    // Forward (X) and backward (P) propagation matrices, as a function of 
    // time. Each matrix is unitary of the form [a -b*; b a*], and can
    // therefore be represented by just the first column, [a; b], which
    // is done here.
    cx_mat X(2, N+1); 
    cx_mat P(2, N+1); 
    cx_mat C(2,2);
    cx_mat Ix(2,2);
    cx_mat Iy(2,2);
    
    Ix(0,0) = complex<double>(0,0); Ix(0,1) = complex<double>(1,0); Ix(1,0) = complex<double>(1,0); Ix(1,1) = complex<double>(0,0);
    Iy(0,0) = complex<double>(0,0); Iy(0,1) = complex<double>(0,-1); Iy(1,0) = complex<double>(0,1); Iy(1,1) = complex<double>(0,0);
    
    // Create output RF
    plhs[0] = mxCreateDoubleMatrix(1, N, mxCOMPLEX);
    RFX = mxGetPr(plhs[0]);
    RFY = mxGetPi(plhs[0]);
    if (!mxIsComplex(prhs[0])) {
        for (t=0; t<N; t++) {RFY[t] = 0; RFX[t] = RFX0[t];}
    } else {
        for (t=0; t<N; t++) {RFY[t] = RFY0[t]; RFX[t] = RFX0[t];}
    }
    
    // Calculate quantities of interest
    dt = tp/N;
    curCost = 0;
    prevCost = 0;
    numOffsets = mxGetNumberOfElements(prhs[2]);
    numB1s = mxGetNumberOfElements(prhs[3]);
    
    // Optimal control main loop
    while (!isConverged)
    {
        idxIter++;
        if (isDebug)  mexPrintf("Iteration %d \n", idxIter);
        for (b = 0; b<numB1s; b++)
        {
            curCost = 0;
            for (f = 0; f<numOffsets; f++)
            {
                if (isDebug) mexPrintf("   f = %d, b = %d \n", f, b);
                // Retrieve target propagator. The target propagator is stored in a cell array of the form
                // U{b, f} = 2x1 complex vector = [a; b]
                // Equivalent to the unitary matrix [a -b*; b a*]
                subs0[0] = b; subs0[1] = f;
                ndims = 2;
                if (mxIsComplex(mxGetCell(prhs[4], mxCalcSingleSubscript(prhs[4], ndims, subs0))))
                {
                    UTRe = mxGetPr(mxGetCell(prhs[4], mxCalcSingleSubscript(prhs[4], ndims, subs0)));
                    UTIm = mxGetPi(mxGetCell(prhs[4], mxCalcSingleSubscript(prhs[4], ndims, subs0)));
                } else {
                    UTRe = mxGetPr(mxGetCell(prhs[4], mxCalcSingleSubscript(prhs[4], ndims, subs0)));
                    UTIm = twoZeros;
                }
                targetU(0) = complex<double>(UTRe[0], UTIm[0]);
                targetU(1) = complex<double>(UTRe[1], UTIm[1]);
                
                // Forward propagation: calculate X(j) and store results.
                X(0,0) = complex<double>(1, 0);
                X(1,0) = complex<double>(0, 0);
                if (isDebug) mexPrintf("   About to forward propagate ... \n");
                for (t = 0; t<N; t++)
                {
                    // Calculate P(1) by forward propagation (a bit slower, but saves memory)
                    Bx = RFX[t]*B1Vec[b]; // kHz
                    By = RFY[t]*B1Vec[b]; // kHz
                    Bz = offsetVec[f]; // kHz
                    Beff = sqrt(Bx*Bx+By*By+Bz*Bz); // kHz
                    // Calculate rotation axis and angle
                    if (Beff==0)
                    { 
                        nx = 0; ny = 0; nz = 1; phi = 0; 
                    } else {
                        nx = Bx/Beff;
                        ny = By/Beff;
                        nz = Bz/Beff;
                        phi = -2*pi*Beff*dt; // LH Rotation
                    }
                    cp = cos(phi/2); sp = sin(phi/2);
                    // Calculate and store X(j) = U(j)*X(j-1) as vector [a; b]
                    U(0,0) = complex<double>(cp, -nz*sp);
                    U(1,0) = complex<double>(ny*sp, -nx*sp);
                    U(0,1) = complex<double>(-ny*sp, -nx*sp);
                    U(0,0) = complex<double>(cp, nz*sp);
                    X.col(t+1) = U*X.col(t);
                    if (isDebug) mexPrintf("      t = %d, B=[%.3f, %.3f, %.3f],  Beff = %.4f, phi = %.4f, n=[%.2f, %.2f, %.2f]\n", t,  Bx, By, Bz, Beff, phi, nx, ny, nz);
                }
                if (isDebug) mexPrintf("   Preparing for backward propagation \n");
                // Backward propagation: calculate P(j) and update field
                // Calculate P(0) = XN!*UF
                P.col(0) = htrans(X.col(N))*targetU;
                if (isDebug) mexPrintf("   Entering backward propagation \n");
                for (t = 0; t<N; t++)
                {
                    Bx = RFX[t]*B1Vec[b]; // kHz
                    By = RFY[t]*B1Vec[b]; // kHz
                    Bz = offsetVec[f]; // kHz
                    Beff = sqrt(Bx*Bx + By*By + Bz*Bz); // kHz
                    // Calculate rotation axis and angle
                    if (Beff==0)
                    {
                        nx = 0; ny = 0; nz = 1; phi = 0;
                    } else {
                        nx = Bx/Beff;
                        ny = By/Beff;
                        nz = Bz/Beff;
                        phi = -2*pi*Beff*dt; // LH Rotation
                    }
                    cp = cos(phi/2); sp = sin(phi/2);
                    U(0,0) = complex<double>(cp, -nz*sp);
                    U(1,0) = complex<double>(ny*sp, -nx*sp);
                    U(0,1) = complex<double>(-ny*sp, -nx*sp);
                    U(0,0) = complex<double>(cp, nz*sp);
                    P.col(t+1) = U*P.col(t);
                    // Calculate the correction to the field
                    fullX(0,0) = X(0, t+1);
                    fullX(1,0) = X(1, t+1);
                    fullX(0,1) = -conj(X(1, t+1));
                    fullX(1,1) = conj(X(0, t+1));
                    
                    fullP(0,0) = P(0, t+1);
                    fullP(1,0) = P(1, t+1);
                    fullP(0,1) = -conj(P(1, t+1));
                    fullP(1,1) = conj(P(0, t+1));
                    
                    // c1 = htrans(X.col(t+1))*P.col(t+1) + htrans(P.col(t+1))*X.col(t+1);
//                     c1 = real(trace(htrans(fullX)*fullP));
//                     c2 = real(trace(htrans(fullP)*Ix*fullX));
//                     c3 = real(trace(htrans(fullP)*Iy*fullX));
//                     RFX[t] = RFX[t] - epsilon*dt*(c2*c1);
//                     RFY[t] = RFY[t] - epsilon*dt*(c3*c1);

                    // Cost
                    curCost = curCost + abs(c1)*abs(c1)/2/N;
                }
                        
            }
        }
        for (t = 0; t<N; t++)
        {
            magB1 = sqrt(RFX[t]*RFX[t] + RFY[t]*RFY[t]);
            if (magB1>maxB1)
            {
                RFX[t] = RFX[t]/magB1*maxB1;
                RFY[t] = RFY[t]/magB1*maxB1;
            }
        }
        curCost = curCost/numOffsets/numB1s/2;
        if (isOutput) mexPrintf("Iteration %d, cost %.5f \n", idxIter, curCost);
        curCost = 0;
        if (idxIter>=maxIter) isConverged=1;
    }
    
}