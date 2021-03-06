function pulse = PulseCreate_eSNOB(center, tp, Nt);
% SYNTAX: pulse = PulseCreate_eSNOB(center, tp, Nt);
%
% Inputs:
%   center - Center of excitation, in kHz
%   tp     - Pulse duration, in ms
%   Nt     - Number of steps
%
% NOTE: Recall Q factor = 2.2 for eSNOB, so FWHM ~ Q/tp
% See also Kupce et. al., JMR series B 106, 300-303 (1995)

Fourier_Coefficients = [0.7500 -0.6176 -0.0373 -0.0005 -0.0182 -0.0058 -0.0036 -0.0051 -0.0031 -0.0017;
                        0      -0.4855  0.1260 -0.0191 -0.0005 -0.0003  0.0017 -0.0013  0.0001 -0.0025];
w = 2*pi/tp;
tt = linspace(0,tp,Nt);
RFamp = 0.*tt;  RFphase = 0.*tt;

for k=1:10
    RFamp = RFamp + 1/tp*(Fourier_Coefficients(1,k)*cos((k-1)*w*tt) + Fourier_Coefficients(2,k)*sin((k-1)*w*tt));
end;
RFphase = 2*pi*center*tt;

pulse.tp = tp;
pulse.RFamp = RFamp;
pulse.RFphase = RFphase;
pulse.Gx = zeros(1,Nt);
pulse.Gy = zeros(1,Nt);
pulse.Gz = zeros(1,Nt);