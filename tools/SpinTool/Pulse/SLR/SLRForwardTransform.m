function [A,B] = SLRForwardTransform(pulse)
% [A,B] = SLRForwardTransform(pulse)
%
% Maps an RF pulse into two polynomials, A and B, as detailed in Pauly's
% 1991 IEEE paper. The full SLR polynomials are then
%
% Ap(z) = A(1) + A(2)*(1/z) + A(3)*(1/z^2) + ... + A(N)*(1/z^(N-1))
% Bp(z) = B(1) + B(2)*(1/z) + B(3)*(1/z^2) + ... + B(N)*(1/z^(N-1))
%
% where N is the number of steps in the pulse, and A(i), B(i) are the
% i-th elements in the vectors A, B, returned by this function. Here
% z = exp(i*w*dt), where w is the offset and dt is the dwell time of the
% pulse.

numSteps = length(pulse.RFamp);
dwellTime = pulse.tp/numSteps;


% The rotation angle applied by the pulse (assuming the hard-pulse 
% approximation, i.e., neglecting the offset.)
phi = zeros(1, numSteps);
phi(1)  = 2*pi*dwellTime*pulse.RFamp(1);

% The Ck and Sk coefficients which appear in Pauly's paper
C1    = cos(phi(1)/2);
S1    = 1i*sin(phi(1)/2)*exp(1i*pulse.RFphase(1));

% These are A1, B1
A = [C1];
B = [S1];

for k=2:numSteps
    phi(k)  = 2*pi*dwellTime*pulse.RFamp(k);
    Ck    = cos(phi(k)/2);
    Sk    = 1i*sin(phi(k)/2)*exp(1i*pulse.RFphase(k));
    Atemp = Ck*[A 0] - conj(Sk)*[0 B];
    Btemp = Sk*[A 0] +       Ck*[0 B];
    A = Atemp;
    B = Btemp;
end