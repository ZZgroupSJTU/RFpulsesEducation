function [fid, spins] = AcquireJ(spins, acqTime, numPoints, Gx, Gy, Gz)
% SYNTAX: [fid, spins] = AcquireJ(spins, acqTime, numPoints, Gx, Gy, Gz)

pulse.tp = acqTime;
pulse.RFamp = zeros(1,numPoints);
pulse.RFphase = zeros(1,numPoints);
pulse.Gx = Gx.*ones(1,numPoints);
pulse.Gy = Gy.*ones(1,numPoints);
pulse.Gz = Gz.*ones(1,numPoints);
isAcquire = 1;

[fid, spins] = PropagateJ(spins, pulse, isAcquire);

% dt = acqTime/numPoints;
% timeVec = [0:dt:(numPoints-1)*dt];
% % Calculate the effective T2 based on the linewidth
% T2 = (1/spins.linewidth/pi)*1000; % ms
% fid = fid.*exp(-timeVec/T2);