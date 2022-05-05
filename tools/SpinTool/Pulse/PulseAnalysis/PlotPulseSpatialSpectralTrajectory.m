function PlotPulseSpatialSpectralTrajectory(pulse, curAxis)
% SYNTAX:
%
%    PlotPulseSpatialSpectralTrajectory(pulse)
%
% Plots the spatial-spectral trajectory and pulse density in the k-t
% plane, using the gradient specified in curAxis, which is a character set
% to 'x', 'y' or 'z' (if omitted, 'z' will be used by default).

if nargin<2
    curAxis = 'z';
end

N = numel(pulse.RFamp);
t = linspace(0, pulse.tp, N);
switch lower(curAxis)
    case 'x'
        k = cumsum(pulse.Gx)*1000; % m^(-1)
    case 'y'
        k = cumsum(pulse.Gy)*1000; % m^(-1)
    case 'z'
        k = cumsum(pulse.Gz)*1000; % m^(-1)
end

dkdt = diff2(k, t);


% The weighting function W is defined in Pauly's 1989 JMR paper "A linear
% class of large tip angle selective excitation pulses" (Eq. 17)
W = pulse.RFamp./abs(dkdt);

figure
subplot(1,2,1)
plot(t, k);
xlabel('ms');
ylabel('m^{-1}');

subplot(1,2,2)
PlotParametricFunction(t,k,W);
zlabel('kHz');
xlabel('ms');
ylabel('m^{-1}');

