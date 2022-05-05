function pulse = PulseCreateFromSweep(duration, amplitude, frequency, gradient)
% SYNTAX:
%
%   pulse = PulseCreateFromSweep(duration, amplitude, frequency, [gradient])
%
% Creates a pulse structure given a duration, amplitude, and instantaneous
% frequency (1/(2*pi)*dphi/dt) of the RF shape.
%
% Input Variables
% Variable Name    Units    Description
% duration         ms       Pulse duration
% amplitude        kHz      RF amplitude
% frequency        kHz      Instantaneous RF frequency = (1/2pi)*(dphi/dt).
%                           If has a different number of elements than 
%                           amplitude, it will be interpolated.
% [gradient]       kHz/mm   An optional 1x3 vector containing the constant
%                           gradient strengths along each of the axes.
%                           Default: [0 0 0]

if nargin<4
    gradient = [0 0 0];
end

N = numel(amplitude);
dt = duration/N;
tt = [0:dt:(N-1)*dt];

% Interpolate frequency if not the same number of elements as amplitude
if numel(frequency)~=N
    Nf = numel(frequency);
    frequency = interp1(linspace(0,1,Nf), frequency, linspace(0,1,N));
end

% Calculate instantaneous phase
ph = cumsum(frequency*2*pi*dt);

pulse.tp = duration;
pulse.RFamp = amplitude;
pulse.RFphase = ph;

pulse.Gx = gradient(1)*ones(1,N);
pulse.Gy = gradient(2)*ones(1,N);
pulse.Gz = gradient(3)*ones(1,N);


