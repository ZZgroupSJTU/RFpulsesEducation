#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Syntax: spinsOut = ApplyZRotation3D(spinsIn, phi)
//
// Description: applies a given LH z-rotationto a 3D spin ensemble, one spin at a time.
//
// phi is the (uniform) rotation angle about the z-axis, in degrees.
// This function mimics ApplyZRotation, but operates on a 3D spin structure
// as initialized using InitSpins3D.
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
    long   idx;
    double *Mo, *Mi, *MoIm, *MiIm;
    long Nx, Ny, Nz;
    
    // Retrieve Pulse data
    phi = *mxGetPr(prhs[1]);  // Deg.
    ct = cos(phi/180.0*pi);
    st = -sin(phi/180.0*pi); // The minus sign makes it a LH rotation
    
    // Retrieve number of spins
    Nx = mxGetN(mxGetField(prhs[0],0,"xVec"));
    Ny = mxGetN(mxGetField(prhs[0],0,"yVec"));
    Nz = mxGetN(mxGetField(prhs[0],0,"zVec"));
    numSpins = Nx*Ny*Nz;
    
    // Create output structure (identical to input structure)
    plhs[0] = mxDuplicateArray(prhs[0]);
    Mi = mxGetPr(mxGetField(prhs[0], 0, "M"));  
    Mo = mxGetPr(mxGetField(plhs[0], 0, "M"));  

    for (idx=0; idx<numSpins; idx++)
    {
        Mo[idx]          = Mi[idx]*ct - Mi[idx+numSpins]*st;
        Mo[idx+numSpins] = Mi[idx]*st + Mi[idx+numSpins]*ct;
    }

    if (mxIsComplex(mxGetField(prhs[0], 0, "M"))==1) {
        MiIm = mxGetPi(mxGetField(plhs[0], 0, "M"));
        MoIm = mxGetPi(mxGetField(plhs[0], 0, "M"));
        for (idx=0; idx<numSpins; idx++)
        {
            MoIm[idx]          = MiIm[idx]*ct - MiIm[idx+numSpins]*st;
            MoIm[idx+numSpins] = MiIm[idx]*st + MiIm[idx+numSpins]*ct;
        }
    };
    
}