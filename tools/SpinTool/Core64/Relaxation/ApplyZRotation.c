#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: spinsOut = ApplyZRotation(spinsIn, phi)
//
// Description: applies a given LH z-rotationto a spin ensemble, one spin at a time.
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
//      phi is the (uniform) rotation angle about the z-axis, in degrees.
//
// Output:
//
//      A spin structure at the end of the pulse.
{
    const double pi = 3.1415926535897932384626433832795;
    int numSpins;  
    int curSpinIdx;
    double vxR, vyR, vzR, vxI, vyI, vzI;      // temporary vector components
    double *Mout, *Min, *MoutIm, *MinIm;
    double phi, ct, st;    // Internal variables for bloch equations
    int    isComplex;                                                     // =1 if magnetization is complex. There are no boolean types in most C standards. Let's use int.
    
    // Retrieve Pulse data
    phi = *mxGetPr(prhs[1]);  // Deg.
    ct = cos(phi/180.0*pi);
    st = -sin(phi/180.0*pi); // The minus sign makes it a LH rotation
    
    // Retrieve number of spins
    numSpins = mxGetNumberOfElements(prhs[0]);
    
    // Create output structure (identical to input structure)
    plhs[0] = mxDuplicateArray(prhs[0]);
    
    for (curSpinIdx=0; curSpinIdx<numSpins; curSpinIdx++)          // Loop over spins
    {
        Mout = mxGetPr(mxGetField(plhs[0], curSpinIdx, "M"));      // Get pointer to magnetization of spin object
        vxR  = Mout[0]; vyR = Mout[1];
        Mout[0] = vxR*ct - vyR*st;
        Mout[1] = vyR*ct + vxR*st;
        // If there is a complex component, do the exact same thing for it
        if (mxIsComplex(mxGetField(prhs[0], curSpinIdx, "M"))==1) {
            MoutIm = mxGetPi(mxGetField(plhs[0], curSpinIdx, "M")); 
            vxI = MoutIm[0]; vyI = MoutIm[1];
            MoutIm[0] = vxI*ct - vyI*st;
            MoutIm[1] = vyI*ct + vxI*st;
        }
    }
    
}