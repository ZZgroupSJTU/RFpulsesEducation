function [spins, fid] = Acquire3D(spins, pulse)
% Acquires an FID from a 3D spin structure. 
%   [spinsOut, fid] = Acquire3D(spins, pulse) Propagates the spins forward
%   while acquiring an FID from the structure. Returns the FID in a complex
%   vector, fid. Note that the RF elements of the input pulse are ignored!

cs = spins.cs;
Gx = pulse.Gx;
Gy = pulse.Gy;
Gz = pulse.Gz;
[xx,yy,zz] = ndgrid(spins.xVec, spins.yVec, spins.zVec);
Mxy = spins.M(:,:,:,1) + 1i*spins.M(:,:,:,2);
Mz = spins.M(:,:,:,3);
numSteps = numel(pulse.RFamp);
dwellTime = pulse.tp/numSteps;
R2 = exp(-dwellTime/spins.T2);

fid = zeros(1, numSteps);

if IsPulseConstGrad(pulse)
    % Constant gradients. This means we can generate the transverse 
    % propgator only once and reuse it at each step.
    Pxy = exp(-1i*2*pi*(xx*Gx(1) + yy*Gy(1) + zz*Gz(1) + cs)*dwellTime)*R2;
    for idxStep=1:numSteps
        fid(idxStep) = sum(Mxy(:));
        Mxy = Mxy.*Pxy;
    end
    Mz = Mz.*exp(-pulse.tp/spins.T1) + (1-exp(-pulse.tp/spins.T1))*spins.M0;
else
    % Non constant gradients. This means we need to re-generate the
    % transverse propagator at every time step. 
    for idxStep=1:numSteps
        
        fid(idxStep) = sum(Mxy);
        Pxy = exp(-1i*2*pi*(xx*Gx(idxStep) + yy*Gy(idxStep) + zz*Gz(idxStep) + cs)*dwellTime)*R2;
        Mxy = Mxy.*Pxy;
    end
    Mz = Mz.*exp(-pulse.tp/spins.T1) + (1-exp(-pulse.tp/spins.T1))*spins.M0;
end

spins.M(:,:,:,1) = real(Mxy);
spins.M(:,:,:,2) = imag(Mxy);
spins.M(:,:,:,3) = Mz;