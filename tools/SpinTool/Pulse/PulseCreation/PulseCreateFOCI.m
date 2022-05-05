function [pulseFOCI, A] = PulseCreateFOCI(pulse, A, G)
% SYNTAX: pulse = PulseCreateFOCI(pulse, envelope, G)
%
% Description: given an adiabatic pulse (e.g. hyperbolic secant), this
% modules the (constant) gradient - assumed along the z-axis - in such
% a way as to increase the BW without increasing the peak RF power 
% (however, SAR will be increased). G is the input (constant) gradient,
% in kHz/mm.
%
% FOCI pulses are described in: 
% Ordidge et al., Magn. Reson. Med. 36:562-6 (1996).
% 
% The envelope input variable, A, is a function assuming values between
% [1, inf], and is equal to 1 in the center of the pulse. It modulates
% B1, the offset, and the gradient as a function of time. The user can
% input one of two envelope types:
% 1. If a number is entered between (0,1), it will determine what
%    fraction of the final B1 shape will be "flat" (see Fig. 3 in Ordidge's
%    paper). The higher the fraction, the wider the BW and the higher
%    the SAR. Note that higher percentages will also mean a larger
%    envelope which might lead to experimentally hard-to-realize shapes
%    for the gradients and/or RF shape.
% 2. Alternatively, a vector having the same number of points as 
%    the RF waveform can be entered. This will be the explicit form of
%    the envelope.

pulseFOCI = pulse;

N = numel(pulse.RFamp);
T = pulse.tp;
dt = T/N;

% Time axis
tt = [0:dt:(N-1)*dt];

% RF amplitude as a function of time
B1 = pulse.RFamp;

% RF offset as a function of time
w = diff(pulse.RFphase)/dt;
w(end+1) = w(end) + (w(end)-w(end-1))/dt*dt;

% Gradient, in mT/m, as a function of time - initially constant
G = ones(1,N)*G;

% Compute shape of A(t), if a percentage is entered
if (numel(A)==1)
    if (A<0.01)
        A = 0.01;
    end
    if (A>1)
        A = 0.5;
    end
    NF = ceil(N*A);
    if (mod(NF,2)==1)
        NF = NF + 1;
    end
    if (NF>N)
        NF = N;
    end
    B1FOCI = B1;
    NC = round(N/2);
    B1FOCI(NC-NF/2:NC+NF/2) = 1; 
    A = B1FOCI./B1;
end

% Adjust pulse amplitude
pulseFOCI.RFamp = B1.*A;

% Adjust pulse phase
wFOCI = w.*A;
phaseFOCI = cumsum(wFOCI)*dt;
phaseFOCI = phaseFOCI - pulse.RFphase(1);
pulseFOCI.RFphase = phaseFOCI;

% Adjust pulse gradient (Store in z-gradient)
pulseFOCI.Gz = G.*A; 

