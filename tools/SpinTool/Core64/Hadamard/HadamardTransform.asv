#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Applies a Hadamard transform to an N-dimensional input array along the
// specified dimension. Currently only defined for 1, 2, 4, 8 and 16 elements.
//
// Syntax: array = HadamardTransform(array, dimension)
//
// Input:
//
//      spins structure:
//      Field           Type                
//      array           n1xn2x...xnM complex
//      dimension       The dimension along which to perform the HT
//
// Output:
//
//      An array of equal dimensions to the input array, HT-ed along the
//      specified dimension. 
//
// Example: suppose a 2*4 matrix has been input
//
//      a00 a01
//      a10 a11
//      a20 a21
//      a30 a31
//
// The different indices correspond to a linear index as follows:
//
//      0   1
//      2   3
//      4   5
//      6   7
//
// If performing a HT over the 1st dimension (the rows), one would have
// to go over sets of rows and HT them
//
{
    int numDimensions;         // Number of dimensions in input array
    int *numElements;          // Number of elements along each dimension
    int idxDim;
    
    // Create output structure (identical to input structure)
    plhs[0] = mxDuplicateArray(prhs[0]);

    numDimensions = mxGetNumberOfDimensions(prhs[0]);
    numElements = mxGetDimensions(prhs[0]);

    for (idxDim=0; idxDim<numDimensions; idxDim++) {
        mexPrintf("numElements[%d] = %d \n", idxDim, numElements[idxDim]);
    }
     

}