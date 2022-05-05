function [fid, spinsOut] = PropagateJSecular(spinsJ, pulse)
% The current function applies a pulse to the given spinsJ, assuming 
% homonuclear case. The spinsJ are a structure of the form:
%
%   spinsJ.rho    = Density matrix of the spinsJ
%   spinsJ.cs     = spinsJ' chemical shifts, in kHz
%   spinsJ.J      = An NxN symmetric matrix containing J-coupling, in Hz.
%   spinsJ.r      = A column vector (matrix) containing the position of the current spin
%   spinsJ.B1     = B1 scaling factor (set to 1 for no B1 inhomogeneity)
%   spinsJ.B0     = B0 inhomogeneity, in kHz (set to 0 for bo B0 inhomogeneity)
%
% The input pulse has a structure:
%
%   pulse.tp      = time of pulse, in ms (#)
%   pulse.RFamp   = amplitude of pulse, in kHz (vector)
%   pulse.RFphase = phase of pulse, in radians (vector)
%   pulse.Gx      \
%   pulse.Gy       |-> Gradients, in kHz/mm (vector)
%   pulse.Gz      /
%
% The current function does not take relaxation into account.

if (size(spinsJ.J,1)~=size(spinsJ.J,2))
    fprintf('PropgateJ error: spinsJ.J is not a square matrix. Aborting!\n');
    beep
    return
end
    

if (length(spinsJ.cs)~=size(spinsJ.J, 1))
    fprintf('PropgateJ error: spinsJ.cs has more elements than spinsJ.J. Aborting!\n');
    beep
    return
end
    
numSpins = length(spinsJ.cs);
Nt = length(pulse.RFamp);
dt = pulse.tp/Nt;

% Step I: Prepare operators
Ix_tot = zeros(2^numSpins, 2^numSpins);
Iy_tot = zeros(2^numSpins, 2^numSpins);
Iz_tot = zeros(2^numSpins, 2^numSpins);
for k=1:numSpins
    Ix_tot = Ix_tot + IxN(k,numSpins);
    Iy_tot = Iy_tot + IyN(k,numSpins);
    Iz_tot = Iz_tot + IzN(k,numSpins);
end;
Ixy_tot = Ix_tot + 1i*Iy_tot;

% Step II: Compute time independent Hamiltonian, H = (chemical shifts) + (J-couplings), in kHz*rad. Recall J is in Hz, not kHz. hbar = 1.
Hcs = zeros(2^numSpins, 2^numSpins);
HJ  = zeros(2^numSpins, 2^numSpins);
for p=1:numSpins
    Hcs = Hcs + 2*pi*spinsJ.cs(p)*IzN(p,numSpins);
end;
% Add J coupling. Use only the upper right triangular half of the matrix spinsJ.J
for p=1:numSpins
    for k=p+1:numSpins
        HJ= HJ + 2*pi*spinsJ.J(p,k)*0.001*IzN(p,numSpins)*IzN(k,numSpins);
    end;
end;

HSpins = Hcs + HJ;

% Step III: Compute time dependent fields
% Compute the operator that will be added to the Hamiltonian for each step in the pulse for the RF (grd. not included!)
H_RF{Nt} = 0.*Iz_tot;
for idxTime=1:Nt
    H_RF{idxTime} = 2*pi*pulse.RFamp(idxTime)*cos(pulse.RFphase(idxTime))*Ix_tot + 2*pi*pulse.RFamp(idxTime)*sin(pulse.RFphase(idxTime))*Iy_tot;
end;

% Step III: Propagate
fid = zeros(1,Nt);
for idxTime=1:Nt
    % Acquire FID
    fid(idxTime) = fid(idxTime) + trace(spinsJ.rho*Ixy_tot);
    % Compute gradient fields: 2*pi*gamma*[Gx(t)*x + Gy(t)*y + Gz(t)*z]
    grdField = 2*pi*(pulse.Gx(idxTime)*spinsJ.r(1) ...
                   + pulse.Gy(idxTime)*spinsJ.r(2) ...
                   + pulse.Gz(idxTime)*spinsJ.r(3))*Iz_tot;
    U = expm(-1i*(HSpins + H_RF{idxTime} + grdField)*dt);
    spinsJ.rho = U*spinsJ.rho*U';
end;
spinsOut = spinsJ;
