function [spinsOut, fid] = AcquireGradRelax(spinsIn,acqTime,numAcqPoints,Gx,Gy,Gz)
% Acquires the fid in the presence of a constant gradients Gx, Gy, Gz, given in kHz/mm.
% Relaxation is taken into account.
%
% Input Parameters:
%   Name            Type         Units            Description      
%   spinsIn         1D spins     N/A              InumAcqPointsut spins
%   acqTime         double       milliseconds     Acquisition duration [1]
%   numAcqPoints    double       Unitless         Number of points
%   Gx              double       kHz/mm           Acquisition gradient (can also be negative)
%   Gy              double       kHz/mm           Acquisition gradient (can also be negative)
%   Gz              double       kHz/mm           Acquisition gradient (can also be negative)

pulse.tp      = acqTime;
pulse.RFamp   = zeros(1,numAcqPoints);
pulse.RFphase = zeros(1,numAcqPoints);
pulse.Gx      = Gx.*ones(1,numAcqPoints);
pulse.Gy      = Gy.*ones(1,numAcqPoints);
pulse.Gz      = Gz.*ones(1,numAcqPoints);

[spinsOut, fid] = AcquirePulseRelax(spinsIn, pulse);