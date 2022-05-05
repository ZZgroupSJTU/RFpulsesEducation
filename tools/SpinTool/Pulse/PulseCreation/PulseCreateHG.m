function pulse = PulseCreateHG(duration, maxAmp, threshold, numSteps, a, b, c, d) 
% Creates a hypergeometric pulse. At high enough B1s it serves as an
% adiabatic pulse. Based on
%
%   Rosenfeld et. al., Analytic solutions of the Bloch equation involving 
%   asymmetric amplitude and frequency modulations. Phys. Rev. A, 54(3):
%   2439-2443, 1996.
%
% Input Variables
% Variable Name     Units     Description
% duration          ms        Total pulse duration
% maxAmp            kHz       Maximal pulse amplitude
% threshold         0-1       The cutoff value, with respect to the max.,
%                             used to determine the "amount" of pulse used
% numSteps          -         
% a,b,c,d                     Pulse parameters:
%                             a>-b the asymmetry parameter
%                             b>0, c>0, d>0 rate parameters
% 
%
% Output Variables
% Variable Name     Units     Description
% pulse             -         The output pulse


% Find the maximal value of R(z)=RMax and the thresholded R(z)=RMin,
% and consequently the range of values of z (and consequently of time)
% Here R(z) is taken from Eq. (10). z is calculated by solving R(z)=RMin,
% and t(z) is given by Eq. (11)
RMax = 1/(2*sqrt(b*(a+b)));
RMin = threshold*RMax;
%zi = (1-2*a*b*RMin^2-sqrt(1-4*a*b*RMin^2 - 4*b^2*RMin^2))/(2*(1+a^2*RMin^2));
%zf = (1-2*a*b*RMin^2+sqrt(1-4*a*b*RMin^2 - 4*b^2*RMin^2))/(2*(1+a^2*RMin^2));
%ti = log(zi^b/(1-zi)^(b+a));
%tf = log(zf^b/(1-zf)^(b+a));
ti = - duration/2;
tf = duration/2;

% Create "unitless" time axis
timeAxis = linspace(ti, tf, numSteps);

RFamp = zeros(1,numSteps);
RFphase = zeros(1,numSteps);
options = optimset('TolX', 1e-45);
for idx=1:numSteps
    z = t2z(timeAxis(idx),a,b, options);
    RFamp(idx) = sqrt(z*(1-z))/(a*z+b);
    RFphase(idx) = log(z^d/(1-z)^(c+d));
end

pulse.tp = duration;
pulse.RFamp = RFamp.*maxAmp;
pulse.RFphase = RFphase;
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);




% Converts between the time variable (t) and the nonlinear z-variable
% defined in Rosenfeld's paper (Eq. 11), ranging from 0 to 1. a and b
% are the same parameters a, b defined in the main function.
function z=t2z(t, a, b, options)
z = fminbnd(@(z) abs(log(z^b/(1-z)^(a+b))-t), 0, 1, options);