function pulse = PulseCreateRosenfeld(a,b,c,d,omega0,numSteps, threshold)
% Creates an asymmetric excitation pulse based on 
%
%   Rosenfeld, Panfil and Zur, Phys. Rev. A: 54(3):2439-43 (1996): Analytic
%   solutions of the Bloch equation involving asymmetric amplitude and
%   frequency modulations
%
% Input Variables
% Name       Description
% a,b,c,d    Variables described in Rosenfeld's paper. In particular,
%            they must uphold:
%            b>0
%            a>-b
% omega0     Scaling coefficient, as appears in above papers. NOT equal to
%            maximum RF power (neither in kHz nor in kHz*rad).
% numSteps   Total number of pulse steps
% threshold  Optional. If not input, set to 0.01. This is the fraction of
%            the maximal RF amplitude at which the RF waveform will be
%            truncated.
%
% Output Variables
% Name       Description
% pulse      The output pulse structure

if nargin<7
    threshold = 0.01;
end

% Find the values of z(t) at the leftmost and rightmost points where the
% RF amplitude reaches (threshold) times its maximal value. Essentially
% what I'm doing here is solving Eq. 10 for 
%   Rmax*threshold = sqrt(z*(1-z))/(a*z+b)
maxAmp = 1/(2*sqrt(b*(a+b))); % unitless (Eq. 10 in Phys Rev A)
cutoffAmp = threshold*maxAmp; 
% Coefficients of quadratic equation for z (note: these are NOT a, b, c)
aa = 1 + cutoffAmp^2*a^2;
bb = 2*a*b*cutoffAmp^2 - 1;
cc = cutoffAmp^2*b^2;
D = sqrt(bb^2 - 4*aa*cc); % Discriminant
zMin = - bb/2 - D/2;
zMax = - bb/2 + D/2;
% Find the corresponding time points (Eq. 11 in Phys Rev A)
tMin = log(zMin^b/(1-zMin)^(a+b));
tMax = log(zMax^b/(1-zMax)^(a+b));
tt = linspace(tMin, tMax, numSteps);

% This takes a few milliseconds but saves the hassle of interpolating, and
% is more exact (assuming the tolerance is set small enough)
zt = zeros(1, numSteps);
options = optimset('TolX', 1e-18, 'MaxIter', 1600);
for idx=1:numSteps
    zt(idx) = fminbnd(@(z) abs(log(z^b./(1-z)^(a+b))-tt(idx)), 0, 1, options);
end

% Calculate the pulse's phase and amplitude
pulse.tp = tMax-tMin;
pulse.RFamp = (omega0/2/pi)*sqrt(zt.*(1-zt))./(a*zt+b);
pulse.RFphase = log((zt.^d)./(1-zt).^(c+d));
pulse.Gx = zeros(1, numSteps);
pulse.Gy = zeros(1, numSteps);
pulse.Gz = zeros(1, numSteps);

