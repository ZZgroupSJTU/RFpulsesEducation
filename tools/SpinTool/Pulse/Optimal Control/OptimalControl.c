#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: RF = OptimalControl(RF0, tp, offsetVec, B1Vec, targetU, maxIter, epsilon, isOutput, maxB1, SAR)
//
//      Variable Name   Type            Units     Description
//
// Input:
//      RF0        1xN complex          kHz       Initial RF guess
//      tp         1x1 double           ms        Pulse duration (fixed)
//      offsetVec  1xM double           kHz       offsets to optimize for
//      B1Vec      1xP double           [scale]   B1 values to optimize for
//      targetU    {PxM}x(2x2) complex  -         Target propagator, as func. of offset & B1
//      maxIter    1x1 int              -         Max. number of iterations
//      epsilon    1x1 double           -         Step size, should be >0 and small
//      isOutput   1x1 boolean (0, 1)   -         Output feedback (if set to 1)
//      maxB1      1x1 double           kHz       Max. RF amplitude
//      SAR        1x1 double           -         SAR penalty constant. Set to
//                                                0 for no penalty, >0 for some penalty.
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
// ! = dagger, UF = target propagator
//               
//   UF <---- UN!*UF ...     <--- U3!*...*UN!*UF <--- U2!*...*UN!*UF <--- U1!*U2!*...*UN!*UF     
//   PN       P(N-1)              P2=U2*P1            P1=U1*P0            P0=XN!*UF
//             =U(N-1)*P(N-2)                                                        
//
// 
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
	double gamma = 42.5774806; // kHz/mT = MHz/T
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
    double SAR = 0.0;
    int subs0[2];
    mwSize ndims = 2;
    
    // Initialize auxiliary structures, of size (N+1)*2
    double **XRe;
    double **XIm;
    double **PRe;
    double **PIm;

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
    SAR = *mxGetPr(prhs[9]);
    
    // Dynamically allocate X, P arrays
    XRe = malloc((N+1)*sizeof(double *));
    XIm = malloc((N+1)*sizeof(double *));
    PRe = malloc((N+1)*sizeof(double *));
    PIm = malloc((N+1)*sizeof(double *));
    for (i=0; i<(N+1); i++)
    {
        XRe[i] = malloc(2*sizeof(double));
        XIm[i] = malloc(2*sizeof(double));
        PRe[i] = malloc(2*sizeof(double));
        PIm[i] = malloc(2*sizeof(double));
    };
    
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
    
    if (isDebug) mexPrintf("numOffsets = %d,  numB1s = %d \n", numOffsets, numB1s);
    
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
                
                // Forward propagation: calculate X(j) and store results.
                XRe[0][0] = 1; XRe[0][1] = 0; XIm[0][0] = 0; XIm[0][1] = 0;  // X(0) = identity
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
                    // XRe[t][0] = Re(a), XIm[t][0] = Im(a)
                    // XRe[t][1] = Re(b), XIm[t][1] = Im(b)
                    XRe[t+1][0] = cp*XRe[t][0] + sp*( nz*XIm[t][0] + nx*XIm[t][1] - ny*XRe[t][1]);
                    XIm[t+1][0] = cp*XIm[t][0] + sp*(-nz*XRe[t][0] - nx*XRe[t][1] - ny*XIm[t][1]); 
                    XRe[t+1][1] = cp*XRe[t][1] + sp*( nx*XIm[t][0] + ny*XRe[t][0] - nz*XIm[t][1]); 
                    XIm[t+1][1] = cp*XIm[t][1] + sp*(-nx*XRe[t][0] + ny*XIm[t][0] + nz*XRe[t][1]);
                    // mexPrintf("      t = %d, B=[%.3f, %.3f, %.3f],  Beff = %.4f, phi = %.4f, n=[%.2f, %.2f, %.2f],  a = %.5f+1i*%.5f,  b = %.5f+1i*%.5f % \n", t,  Bx, By, Bz, Beff, phi, nx, ny, nz, XRe[t+1][0], XIm[t+1][0], XRe[t+1][1], XIm[t+1][1]);
                }
                if (isDebug) mexPrintf("   Cur offset = %.3f, X[N][0] at end of forward propagation: %.5f + i*%.5f \n", offsetVec[f], XRe[N][0], XIm[N][0]);
                if (isDebug) mexPrintf("   Preparing for backward propagation \n");
                // Backward propagation: calculate P(j) and update field
                // Calculate P(0) = XN!*UF
                // XN  = [a -b*]
                //       [b  a*]
                // XN!*UF = [a* b*]*[c -d*] = [a*c+b*d   -  ]
                //          [-b a ] [d  c*] = [-bc+ad    -  ] (here * = conjugation)
                // If C=A*B then (here * = multiplication):
                // ReC = A[0]*B[0] - A[1]*B[1] = ReA*ReB - ImA*ImB
                // ImC = A[0]*B[1] + A[1]*B[0] = ReA*ImB + ImA*ReB
                PRe[0][0] = ( XRe[N][0]*UTRe[0] + XIm[N][0]*UTIm[0]) 
                          + ( XRe[N][1]*UTRe[1] + XIm[N][1]*UTIm[1]);
                PIm[0][0] = ( XRe[N][0]*UTIm[0] - XIm[N][0]*UTRe[0])
                          + ( XRe[N][1]*UTIm[1] - XIm[N][1]*UTRe[1]);
                PRe[0][1] = (-XRe[N][1]*UTRe[0] + XIm[N][1]*UTIm[0]) 
                          + ( XRe[N][0]*UTRe[1] - XIm[N][0]*UTIm[1]);
                PIm[0][1] = (-XRe[N][1]*UTIm[0] - XIm[N][1]*UTRe[0])
                          + ( XRe[N][0]*UTIm[1] + XIm[N][0]*UTRe[1]);
                // mexPrintf("U: a = %.3f + 1i*%.3f,  b = %.3f+1i*%.3f \n", UTRe[0], UTIm[0], UTRe[1], UTIm[1]);
                // mexPrintf("X: a = %.3f + 1i*%.3f,  b = %.3f+1i*%.3f \n", XRe[N][0], XIm[N][0], XRe[N][1], XIm[N][1]);
                // mexPrintf("P: a = %.3f + 1i*%.3f,  b = %.3f+1i*%.3f \n", PRe[0][0], PIm[0][0], PRe[0][1], PIm[0][1]);
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
                    PRe[t+1][0] = cp*PRe[t][0] - sp*(-nz*PIm[t][0] - nx*PIm[t][1] + ny*PRe[t][1]);
                    PIm[t+1][0] = cp*PIm[t][0] - sp*( nz*PRe[t][0] + nx*PRe[t][1] + ny*PIm[t][1]); 
                    PRe[t+1][1] = cp*PRe[t][1] - sp*(-nx*PIm[t][0] - ny*PRe[t][0] + nz*PIm[t][1]); 
                    PIm[t+1][1] = cp*PIm[t][1] - sp*( nx*PRe[t][0] - ny*PIm[t][0] - nz*PRe[t][1]);
                    // Calculate the correction to the field
                    // Let
                    // X = [a -b*]   P = [c -d*]
                    //     [b  a*]       [d  c*]
                    // Then
                    //
                    //            a5*a1            a5*a2            a6*a1            a6*a2
                    // dphi/du1 ~ Im[bc*]Re[ac*] + Im[bc*]Re[bd*] + Im[ad*]Re[ac*] + Im[ad*]Re[bd*]
                    //
                    //            a3*a1            a3*a2            a4*a1            a4*a2
                    // dphi/du2 ~ Re[bc*]Re[ac*] + Re[bc*]Re[bd*] - Re[ad*]Re[ac*] - Re[ad*]Re[bd*]
                    // Thus:
                    // Compute a1=Re[ac*], a2=Re[bd*], a3=Re[bc*], a4=Re[ad*], a5=Im[bc*], a6=Im[ad*]
                    a1 =  XRe[t+1][0]*PRe[t+1][0] + XIm[t+1][0]*PIm[t+1][0];
                    a2 =  XRe[t+1][1]*PRe[t+1][1] + XIm[t+1][1]*PIm[t+1][1];
                    a3 =  XRe[t+1][1]*PRe[t+1][0] + XIm[t+1][1]*PIm[t+1][0];
                    a4 =  XRe[t+1][0]*PRe[t+1][1] + XIm[t+1][0]*PIm[t+1][1];
                    a5 = -XRe[t+1][1]*PIm[t+1][0] + XIm[t+1][1]*PRe[t+1][0];
                    a6 = -XRe[t+1][0]*PIm[t+1][1] + XIm[t+1][0]*PRe[t+1][1];

                    RFX[t] = RFX[t] - epsilon*dt*((a5+a6)*(a1+a2)) - 2*epsilon*dt*SAR*RFX[t];
                    RFY[t] = RFY[t] - epsilon*dt*((a4-a3)*(a1+a2)) - 2*epsilon*dt*SAR*RFY[t];

                    // Cost
                    q1 = 2*(a1+a2);
                    curCost = curCost + fabs(q1)*fabs(q1)/2/N;
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
    
    // Free memory
    for (i=0; i<(N+1); i++)
    {
        free(XRe[i]);
        free(XIm[i]);
        free(PRe[i]);
        free(PIm[i]);
    };
    free(XRe);
    free(XIm);
    free(PRe);
    free(PIm);
    
}