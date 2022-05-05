function outputMatrix = KroneckerProduct(varargin)
% Generates the Kronecker product of multiple entry matrices. 
% (The built-in MATLAB kron function merely operates on two inputs.
% this extends it.)
%
% Inputs:
%
% Any set of matrices
%
% Output:
%
% A matrix which is the kronecker product of the input matrices.
%
% Example:
%
% myOutputMatrix = KroneckerProduct(myInputMatrix1, myInputMatrix2, myInputMatrix3)

numInputMatrices = nargin;
outputMatrix = kron(varargin{numInputMatrices-1}, varargin{numInputMatrices});
for idx=1:(numInputMatrices-2)
    outputMatrix=kron(varargin{numInputMatrices-2},outputMatrix);
end;