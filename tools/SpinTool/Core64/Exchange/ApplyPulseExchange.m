function [spins, fid] = ApplyPulseExchange(spins, pulse, pulseOffset, pulsePhase)
% ApplyPulseExchange  Applies a pulse to an exchanging spin structure.
%   [spins, fid] = ApplyPulseExchange(spins, pulse)  Applies a pulse
%   structure to an exchanging homonuclear proton spin structure (see below). 
%
%   [spins, fid] = ApplyPulseExchange(spins, pulse, pulseOffset) Allows for
%   the specification of a constant pulse offset, in kHz. If omitted,
%   set to 0.
%
%   [spins, fid] = ApplyPulseExchange(spins, pulse, pulseOffset, pulsePhase) 
%   pulsePhase is in degrees, and allows for the introduction of a constant
%   pulse phase. If omitted, set to 0.
%
% The spins are in a structure of the form (N: number of spins per position):
%   spins(i).r      [x; y; z] in mm
%   spins(i).M      3N*1 vector of the form [Mx1; My1; Mz1; Mx2; My2;
%                   Mz2; ... ]
%   spins(i).cs     N*1 vector of chemical shifts of spins, in ppm
%   spins(i).T1     N*1 vector of T1s of spins, ms
%   spins(i).T2     N*1 vector of T2s of spins, ms
%   spins(i).M0     N*1 vector of equilibrium magnetization values
%   spins(i).B1     Scales RF at the position
%   spins(i).RS     Receiver sensitivity scaling at position r
%   spins(i).K      N*N exchange matrix, with exchange constants in Hz.
%                   Note the diagonal in each row needs to equal 
%                   -sum(row, except diagonal)
% The pulse structure is as follows:
%   pulse.tp        Pulse duration in ms
%   pulse.RFamp     1xM vector of pulse amplitude, in kHz
%   pulse.RFphase   1xM vector of pulse phase, in radians
%   pulse.Gx        1xM vector of x-gradient
%   pulse.Gy        1xM vector of x-gradient
%   pulse.Gz        1xM vector of x-gradient
%
% is the same old pulse structure used in ApplyPulseRelax and so
% forth. 

if nargin<3, pulseOffset = 0; end
if nargin<4, pulsePhase = 0; end

numSpins = numel(spins);
numSteps = numel(pulse.RFamp);
dt       = pulse.tp/numSteps; % ms
fid      = zeros(1, numSteps);

RFx = pulse.RFamp.*cos(pulse.RFphase + pulsePhase/180*pi); % kHz
RFy = pulse.RFamp.*sin(pulse.RFphase + pulsePhase/180*pi); % kHz

BxMat = [ 0 0  0;
          0 0  1;
          0 -1 0];

ByMat = [ 0 0 -1;
          0 0  0;
          1 0  0];

BzMat = [ 0  1 0;
          -1 0 0;
           0 0 0];


for idxSpin=1:numSpins
    numSites = numel(spins(idxSpin).cs);
    
    Bx = kron(eye(numSites), BxMat);
    By = kron(eye(numSites), ByMat);
    Bz = kron(eye(numSites), BzMat);

    E = eye(3*numSites);
    
    % Create a vector of relaxation values (T1, T2) for each site, in ms
    if numel(spins(idxSpin).T1)~=numSites 
        if numel(spins(idxSpin).T1)==1
            T1 = repmat(spins(idxSpin).T1, numSites, 1);
        else
            error('Number of T1 values is incompatible.');
        end
    else
        T1 = spins(idxSpin).T1;
    end
    if numel(spins(idxSpin).T2)~=numSites 
        if numel(spins(idxSpin).T2)==1
            T2 = repmat(spins(idxSpin).T2, numSites, 1);
        else
            error('Number of T2 values is incompatible.');
        end
    else
        T2 = spins(idxSpin).T2;
    end
    if numel(spins(idxSpin).M0)~=numSites 
        if numel(spins(idxSpin).M0)==1
            M0 = repmat(spins(idxSpin).M0, numSites, 1);
        else
            error('Number of M0 values is incompatible.');
        end
    else
        M0 = spins(idxSpin).M0;
    end
    if isrow(M0), M0 = M0.'; end 
    
    % Create the transformation matrix to the "second rotating frame"
    % of the pulse, in case it has an offset (in this frame, the offset
    % related time dependence is nulled). 
    offsetAngle = 2*pi*pulseOffset*dt; % in radians
    U           = kron(eye(numSites), RotMat([0 0 1], -offsetAngle));
    effField    = 2*pi*pulseOffset*Bz; % in rad*kHz
    
    % Create 3N*3N relaxation matrix, in kHz
    RDiag = [];
    for idxSite=1:numSites
        RDiag = [RDiag, 1/T2(idxSite), 1/T2(idxSite), 1/T1(idxSite)];
    end
    R = -diag(RDiag);
    
    % Create 3N*1 vector of equilibrium magnetization 
    M0Vec = kron(M0, [0;0;1]);
    
    % Create 3N*3N exchange matrix, in kHz
    K = kron(spins(idxSpin).K, eye(3));
    
    % Create 3N*3N chemical shift matrix, CS, in rad*kHz
    % This includes any effective field brought upon by the constant pulse offset.
    CS = zeros(3*numSites, 3*numSites);
    for idxSite=1:numSites
        CS((idxSite-1)*3+1, (idxSite-1)*3+2) = 2*pi*(spins(idxSpin).cs(idxSite) - pulseOffset);  % kHz
        CS((idxSite-1)*3+2, (idxSite-1)*3+1) = -2*pi*(spins(idxSpin).cs(idxSite) - pulseOffset); % kHz
    end
    
    % Calculate time evolution in the original (non rotating) frame
    for idxTime=1:numSteps
        fid(idxTime) = fid(idxTime) + spins(idxSpin).RS*sum(spins(idxSpin).M(1:3:end)) + 1i*sum(spins(idxSpin).M(2:3:end));
        RF = 2*pi*Bx*RFx(idxTime) + 2*pi*By*RFy(idxTime);
        W  = RF - CS;
        T = R + K + W;
        P = U*expm(T*dt);
        spins(idxSpin).M = P*spins(idxSpin).M + (E - P)*inv(T)*R*M0Vec;
    end
    
end

