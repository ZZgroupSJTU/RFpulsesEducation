function [fid, spins] = PropagateJ(spins, pulse, acqType, affectedNuclei, freqRange, gradUnits, RFUnits)
% SYNTAX: [fid, spinsOut] = PropagateJ(spins, pulse, isAcquire, affectedNuclei, freqRange, gradUnits, RFUnits)
%
% Input Variables
% Variable Name    Description
% spins            Input spin structure. For more details, see below.
% pulse            Input pulse structure. For more details, see below.
% affectedNuclei   A cell array indicating which nuclei and ppm ranges
%                  are to be affected. For example:
%                  affectedNuclei = {{'1h', 3, 4}, {'13c'}}
%                  This indicates all 13c nuclei will be acted on,
%                  as well as any protons with chemical shifts between
%                  3 and 4 ppm (before taking into account spins.csCenter).
%                  This is a "cheap" way of simulating (perfect) selective 
%                  pulses.
% acqType          0: No acquisition.
%                  1: Acquire on, sum FID for all molecules.
%                  2: Acquire on, produce cell array of FIDs for each of the molecules.
% freqRange        Used to crudely simulate frequency selective pulses. 
%                  Any nucleus falling outside the given frequency range
%                  will not be affected by the pulse. freqRange is a 
%                  1x2 vector of [lowest ppm, highest ppm]
% gradUnits        'mt/m' or 'khz/mm'. Default: khz/mm, assuming protons.
%                  This tells the software in which units the Gx, Gy, Gz
%                  fields of the pulse are given.
% RFUnits          'uT' or 'kHz'. Default: kHz, assuming protons. This tells
%                  the software in which units pulse.RFamp is given. 
%
%
% The current function applies a pulse to the given spins, assuming 
% homonuclear case. The spins are a structure of the form specified in
% InitSpinsJ. The input pulse has a structure:
%
%   pulse.tp      = time of pulse, in ms (#)
%   pulse.RFamp   = amplitude of pulse, in microTesla (vector)
%   pulse.RFphase = phase of pulse, in radians (vector)
%   pulse.Gx      \
%   pulse.Gy       |-> Gradient vectors, in gradUnits
%   pulse.Gz      /


% =========================================================================
% Check inputs
% =========================================================================

isLHRot = 1; % If set to 1, the left hand convention for rotations will be used throughout

isSecular = spins.isSecular;
numMolecules = numel(spins.molecule);
if (nargin<3), acqType = 0; end
if (nargin<4)
    affectedNuclei = [];
end
if (nargin<5)
    freqRange = [];
end
if nargin<6
    % Assume by default pulse.Gx, pulse.Gy, pulse.Gz are in kHz/mT FOR PROTONS
    gradUnits = 'khz/mm';
end
if nargin<7
    % Assume by default pulse.RFamp is in kHz FOR PROTONS
    RFUnits = 'khz';
end

switch lower(gradUnits)
    case 'khz/mm'
        % Since this simulation expects gradient values in mT/m, convert
        % to mT/m first
        pulse.Gx = pulse.Gx*1000/GetGyromagneticRatio('1h');
        pulse.Gy = pulse.Gy*1000/GetGyromagneticRatio('1h');
        pulse.Gz = pulse.Gz*1000/GetGyromagneticRatio('1h');
    case 'mt/m'
        % No need to do anything
    otherwise
        error('Unrecognized gradient units %s\n', gradUnits);
end

switch lower(RFUnits)
    case 'khz'
        % Convert from kHz (assuming protons) to uT
        pulse.RFamp = pulse.RFamp*1000/GetGyromagneticRatio('1h'); 
    case 'ut'
        % No need to do anything.
    case 'mt'
        % Convert from mT to uT
        pulse.RFamp = pulse.RFamp*1000; 
    case 't'
        % Convert from T to uT
        pulse.RFamp = pulse.RFamp*1e6;
    otherwise
        error('Unrecognized RF units %s\n', RFUnits);
end

if isempty(affectedNuclei)
    % If affectedNuclei is set to [], all nuclei will be affected by the pulse.
    affectedNuclei = {};
    for idxMolecule=1:numMolecules
        if iscell(spins.molecule(idxMolecule).csVec)
            numSubSystems = numel(spins.molecule(idxMolecule).csVec);
            csVec = spins.molecule(idxMolecule).csVec;
        else
            numSubSystems = 1;
            csVec = {spins.molecule(idxMolecule).csVec};
        end
        for idxSys=1:numSubSystems
            numNuclei = numel(csVec{idxSys});
            affectedNuclei{idxMolecule}{idxSys} = ones(1,numNuclei);
        end
    end
else
    % If affectedNuclei is explicitly specified, let's check the syntax is valid.
    % We will also fix minor inconsistencies if we encounter them.
    if numel(affectedNuclei)~=numMolecules
        error('affectedNuclei has incorrect number of elements, %d (this should equal the number of molecules, which is %d)', numel(affectedNuclei), numMolecules);
    end
    % Check that formatting is consistent; i.e., that each "affectedNuclei"
    % cell is divided into sub-cells describing the various sub-systems
    % (i.e. is "in sync" with the csVec)
    tempAffectedNuclei = [];
    for idxMolecule=1:numMolecules
        % Make sure each sub-system is in its own cell
        if ~iscell(affectedNuclei{idxMolecule})
            tempAffectedNuclei{idxMolecule}{1} = affectedNuclei{idxMolecule};
        else
            tempAffectedNuclei{idxMolecule} = affectedNuclei{idxMolecule};
        end
        % Now check the size of affectedNuclei matches up with the size of csVec
        numAffectedNucleiSystems = numel(tempAffectedNuclei{idxMolecule});
        if iscell(spins.molecule(idxMolecule).csVec)
            numSubSystems = numel(spins.molecule(idxMolecule).csVec);
        else
            numSubSystems = 1;
        end
        if numSubSystems~=numAffectedNucleiSystems
            error('affectedNuclei has a different number of subsystems (%d) compared to cs vector (%d subsystems) for molecule %d (%s)', numAffectedNucleiSystems, numSubSystems, idxMolecule, spins.molecule(idxMolecule).abbrev);
        end
    end 
    affectedNuclei = tempAffectedNuclei;
end

if ~isempty(freqRange)
    for idxMolecule=1:numMolecules
        if iscell(spins.molecule(idxMolecule).csVec)
            % Molecule has multiple subsystems (e.g. NAA CH3 and NAA multiplet). 
            numSubSystems = numel(spins.molecule(idxMolecule).csVec);
            for idxSys=1:numSubSystems
                affectedNuclei{idxMolecule}{idxSys} = affectedNuclei{idxMolecule}{idxSys} & (spins.molecule(idxMolecule).csVec{idxSys}<=freqRange(2)) & (spins.molecule(idxMolecule).csVec{idxSys}>=freqRange(1));       
            end
        else
            % Only one sub-system in the molecule.
            affectedNuclei{idxMolecule}{1} = affectedNuclei{idxMolecule}{1} & (spins.molecule(idxMolecule).csVec<=freqRange(2)) & (spins.molecule(idxMolecule).csVec>=freqRange(1));       
        end
    end
end

Nt = length(pulse.RFamp);
dt = pulse.tp/Nt; % in ms
switch acqType
    case 0
        fid = [];
    case 1
        fid = zeros(1,Nt);
    case 2
        for idx=1:numMolecules
            fid{idx} = zeros(1,Nt);
        end
end

% =========================================================================
% Decide on propagation scheme
% =========================================================================

% If no acquisition is required, and if the gradient of the pulse is
% time independent, and if the spins are "regularly" organized in space
% enough, we can create a 1D simulation (as a function of offset) to save 
% time. Here, we check the conditions for this are met. Basically, we 
% need:
% 1. The gradient is not time dependent and along one axis.
% 2. There is no B0 or B1 inhomogeneity (or, if there are, they are
%    constant across all spins in a molecule).
% 3. This function does not compute an FID, returning an empty set 
%    fid=[] instead.
% 4. Spins are spatially "regular" enough along the gradient's axis.

is1DSim = 1;

idxAxis = 1;
sumGx = sum(abs(pulse.Gx));
sumGy = sum(abs(pulse.Gy));
sumGz = sum(abs(pulse.Gz));
% If no gradients are present, forego the 1D simulation
if (sumGx==0) && (sumGy==0) && (sumGz==0), is1DSim = 0; end 
% If gradients are present along more than one axis, forgo the 1D simulation
if ((sumGx>0)+(sumGy>0)+(sumGz>0))>1, is1DSim = 0; end
% Check along which axis the gradient is applied
if sumGx>0, idxAxis=1; grad=pulse.Gx(1); end
if sumGy>0, idxAxis=2; grad=pulse.Gy(1); end
if sumGz>0, idxAxis=3; grad=pulse.Gz(1); end

B0MoleculeVec = zeros(1, numMolecules);
B1MoleculeVec = zeros(1, numMolecules);

% Check which spatial positions the spins have
for idxMolecule=1:numMolecules
    numSpins = numel(spins.molecule(idxMolecule).spin);
    posVec = zeros(1,numSpins); 
    B0Vec = zeros(1,numSpins);
    B1Vec = zeros(1,numSpins);
    for idxSpin=1:numSpins
        posVec(idxSpin) = spins.molecule(idxMolecule).spin(idxSpin).r(idxAxis); % mm
        B0Vec(idxSpin) =  spins.molecule(idxMolecule).spin(idxSpin).B0;         % ppm
        B1Vec(idxSpin) =  spins.molecule(idxMolecule).spin(idxSpin).B1;         % Scaling
    end
    posVecUnique = unique(posVec);
    B0Vec = unique(B0Vec);
    B1Vec = unique(B1Vec);
    % Only use the simplistic 1D simulation if all spins (for this molecule)
    % have the same B0 and B1 offset.
    if numel(B0Vec)>2, is1DSim = 0; else, B0MoleculeVec(idxMolecule)=B0Vec; end
    if numel(B1Vec)>2, is1DSim = 0; else, B1MoleculeVec(idxMolecule)=B1Vec; end
    rVec{idxMolecule} = posVecUnique; % All unique spin positions for the given molecule
    for idxPos=1:numel(posVecUnique)
        spinPosIndices{idxMolecule, idxPos} = find(posVec==posVecUnique(idxPos));
    end
end

% =========================================================================
%
%
%                                  SIMULATE
%
%
% =========================================================================

if is1DSim
    % =====================================================================
    %
    %                      SIMULATE: 1D PROJECTION
    %
    % =====================================================================
    for idxMolecule=1:numMolecules
        if iscell(spins.molecule(idxMolecule).csVec)
            numSubSystems = numel(spins.molecule(idxMolecule).csVec);
            csVec = spins.molecule(idxMolecule).csVec;
            JMatrix = spins.molecule(idxMolecule).JMatrix;
            nucleus = spins.molecule(idxMolecule).nucleus;
        else
            numSubSystems = 1;
            csVec = {spins.molecule(idxMolecule).csVec};
            JMatrix = {spins.molecule(idxMolecule).JMatrix};
            nucleus = {spins.molecule(idxMolecule).nucleus};
        end
        for idxSys=1:numSubSystems
            numNuclei = numel(nucleus{idxSys});
            numPositions = numel(rVec{idxMolecule});
            [gmRatio, spin] = GetGyromagneticRatio(nucleus{idxSys});
            numDims = prod(spin*2+1);
            % Step I: Prepare operators (weighted by gyromagnetic ratios)
            Ix_tot_RF = zeros(numDims, numDims);
            Iy_tot_RF = zeros(numDims, numDims);
            Iz_tot_grd = zeros(numDims, numDims);
            for k=1:numNuclei
                Ix_tot_RF = Ix_tot_RF + IxN(k,numNuclei,spin)*affectedNuclei{idxMolecule}{idxSys}(k)*gmRatio(k);
                Iy_tot_RF = Iy_tot_RF + IyN(k,numNuclei,spin)*affectedNuclei{idxMolecule}{idxSys}(k)*gmRatio(k);
                Iz_tot_grd = Iz_tot_grd + IzN(k,numNuclei,spin)*gmRatio(k); 
            end

            % ========================================================================
            % Compute CS Hamiltonian (sans spatially dependent B0 offsets)
            % ========================================================================

            Hcs = zeros(numDims, numDims);
            for p=1:numNuclei
                %    <----------------------uT------------------>
                %    <-------------- ppm ------------->  <---T-->       <-kHz/mT->
                cs = (csVec{idxSys}(p) - spins.csCenter)*spins.B0*0.001*gmRatio(p); % in kHz
                Hcs = Hcs - 2*pi*cs*IzN(p,numNuclei,spin); % rad*kHz
            end

            % ========================================================================
            % Compute J-Coupling Hamiltonian, in 2*pi*kHz
            % ========================================================================

            HJ  = zeros(numDims, numDims);
            for p=1:numNuclei
                for k=p+1:numNuclei
                    HJ = HJ + 2*pi*JMatrix{idxSys}(p,k)*0.001*( ...
                         (1-isSecular)*(IxN(p,numNuclei,spin)*IxN(k,numNuclei,spin) + ...
                                        IyN(p,numNuclei,spin)*IyN(k,numNuclei,spin)) ...
                         + IzN(p,numNuclei,spin)*IzN(k,numNuclei,spin)); % rad*kHz
                end
            end
            H0 = Hcs + HJ; % rad*kHz

            % ========================================================================
            % Propagate each molecule through time
            % ========================================================================

            % Calculate the RF Hamiltonian WITH B1 scaling.
            % Remember: pulse.RFamp is in uT. Ix_tot_RF, Iy_tot_RF are in kHz/mT.
            B1Scaling = B1MoleculeVec(idxMolecule)*spins.B1; % (Local B1)*(Global B1)
            for idxTime=1:Nt
                H_RF{idxTime} = 2*pi*pulse.RFamp(idxTime)*0.001*...
                                    (cos(pulse.RFphase(idxTime))*Ix_tot_RF ...
                                   + sin(pulse.RFphase(idxTime))*Iy_tot_RF); % In kHz*rad
            end
            % There is only one offset (by assumption!) in 1D projection simplified simulations
            offset = B0MoleculeVec(idxMolecule)*spins.B0; % uT
            for idxPos=1:numPositions
                % =====================================
                % Compute Time-Dependent RF Hamiltonian 
                % =====================================
                % grad is in mT/m, rVec is in mm, offset is in uT, HGrad is in kHz*rad
                HGrad = 2*pi*(grad*rVec{idxMolecule}(idxPos) + offset)*0.001*Iz_tot_grd; 
                
                U = eye(numDims);
                for idxTime=1:Nt
                    % ==================
                    % Compute propagator
                    % ==================
                    if ((acqType~=0) && isSecular) 
                        P = (-1)^(isLHRot-1)*1i*(H0 + HGrad)*dt; % No RF. dt is in ms. HGrad, H0 is in rad*kHz.
                        U = diag(exp(diag(P)))*U;
                    else
                        U = expm((-1)^(isLHRot-1)*1i*(H0 + HGrad + B1Scaling*H_RF{idxTime})*dt)*U;
                    end
                end

                % =========================================
                % Propagate all spins with a given position
                % =========================================
                numSpins = numel(spinPosIndices{idxMolecule, idxPos});
                for idxSpin=1:numSpins
                    curSpinIdx = spinPosIndices{idxMolecule, idxPos}(idxSpin);
                    if iscell(spins.molecule(idxMolecule).spin(curSpinIdx).rho)
                        spins.molecule(idxMolecule).spin(curSpinIdx).rho{idxSys} = U*spins.molecule(idxMolecule).spin(curSpinIdx).rho{idxSys}*U';
                    else
                        spins.molecule(idxMolecule).spin(curSpinIdx).rho = U*spins.molecule(idxMolecule).spin(curSpinIdx).rho*U';
                    end
                end  
            end  
        end  % Loop over sub-systems
    end   % Loop over molecules
else
    % =====================================================================
    %
    %                         SIMULATE: FULL 3D
    %
    % =====================================================================
    for idxMolecule=1:numMolecules
        if iscell(spins.molecule(idxMolecule).csVec)
            numSubSystems = numel(spins.molecule(idxMolecule).csVec);
            csVec = spins.molecule(idxMolecule).csVec;
            JMatrix = spins.molecule(idxMolecule).JMatrix;
            nucleus = spins.molecule(idxMolecule).nucleus;
        else
            numSubSystems = 1;
            csVec = {spins.molecule(idxMolecule).csVec};
            JMatrix = {spins.molecule(idxMolecule).JMatrix};
            nucleus = {spins.molecule(idxMolecule).nucleus};
        end
        for idxSys=1:numSubSystems
            numNuclei = numel(nucleus{idxSys});
            numSpins = numel(spins.molecule(idxMolecule).spin);
            [gmRatio, spin] = GetGyromagneticRatio(nucleus{idxSys});
            numDims = prod(spin*2+1);
            % Step I: Prepare operators
            Ix_tot = zeros(numDims, numDims);
            Iy_tot = zeros(numDims, numDims);
            Ix_tot_RF = zeros(numDims, numDims);
            Iy_tot_RF = zeros(numDims, numDims);
            Iz_tot_grd = zeros(numDims, numDims);
            for k=1:numNuclei
                Ix_tot = Ix_tot + IxN(k,numNuclei,spin)*gmRatio(k);
                Iy_tot = Iy_tot + IyN(k,numNuclei,spin)*gmRatio(k);
                Ix_tot_RF = Ix_tot_RF + IxN(k,numNuclei,spin)*affectedNuclei{idxMolecule}{idxSys}(k)*gmRatio(k);
                Iy_tot_RF = Iy_tot_RF + IyN(k,numNuclei,spin)*affectedNuclei{idxMolecule}{idxSys}(k)*gmRatio(k);
                Iz_tot_grd = Iz_tot_grd + IzN(k,numNuclei,spin)*gmRatio(k);
            end
            Ixy_tot = Ix_tot + 1i*Iy_tot;

            % ========================================================================
            % Compute CS Hamiltonian (sans spin-dependent B0 offsets), in kHz*rad
            % ========================================================================

            Hcs = zeros(numDims, numDims);
            for p=1:numNuclei
                cs = (csVec{idxSys}(p) - spins.csCenter)*spins.B0*gmRatio(p)/1000; % in kHz
                Hcs = Hcs - 2*pi*cs*IzN(p,numNuclei,spin); % rad*kHz
            end

            % ========================================================================
            % Compute J-Coupling Hamiltonian, in 2*pi*kHz
            % ========================================================================

            HJ  = zeros(numDims, numDims);
            for p=1:numNuclei
                for k=p+1:numNuclei
                    HJ = HJ + 2*pi*JMatrix{idxSys}(p,k)*0.001*( ...
                         (1-isSecular)*(IxN(p,numNuclei,spin)*IxN(k,numNuclei,spin) + ...
                                        IyN(p,numNuclei,spin)*IyN(k,numNuclei,spin)) ...
                         + IzN(p,numNuclei,spin)*IzN(k,numNuclei,spin)); % rad*kHz
                end
            end
            H0 = Hcs + HJ; % rad*kHz
            
            % ========================================================================
            % Propagate each molecule through time
            % ========================================================================

            % Calculate the RF Hamiltonian WITHOUT B1 scaling, in kHz*rad.
            % Remember: pulse.RFamp is in microtesla, while Ix_tot_RF is in kHz/mT.
            for idxTime=1:Nt
                H_RF{idxTime} = 2*pi*pulse.RFamp(idxTime)*0.001*...
                                    (cos(pulse.RFphase(idxTime))*Ix_tot_RF +...
                                     sin(pulse.RFphase(idxTime))*Iy_tot_RF); % In kHz*rad
            end
            for idxSpin=1:numSpins
                % =====================================
                % Compute Time-Dependent RF Hamiltonian
                % =====================================
                B1Scaling = spins.molecule(idxMolecule).spin(idxSpin).B1*spins.B1; % (Local B1)*(Global B1)
                T2 = (1/(spins.linewidth + spins.molecule(idxMolecule).spin(idxSpin).linewidth)/pi)*1000; % Effective T2, in ms (for ad-hoc acquisition line broadening)
                for idxTime=1:Nt
                    curTime = (idxTime-1)*dt; % ms
                    % Acquire FID
                    switch (acqType)
                        case 1
                            if iscell(spins.molecule(idxMolecule).spin(idxSpin).rho)
                                fid(idxTime) = fid(idxTime) + trace(spins.molecule(idxMolecule).spin(idxSpin).rho{idxSys}*Ixy_tot)*exp(-curTime/T2); % Insert ad-hoc line broadening
                            else
                                fid(idxTime) = fid(idxTime) + trace(spins.molecule(idxMolecule).spin(idxSpin).rho*Ixy_tot)*exp(-curTime/T2); % Insert ad-hoc line broadening
                            end
                        case 2
                            if iscell(spins.molecule(idxMolecule).spin(idxSpin).rho)
                                fid{idxMolecule}(idxTime) = fid{idxMolecule}(idxTime) + trace(spins.molecule(idxMolecule).spin(idxSpin).rho{idxSys}*Ixy_tot)*exp(-curTime/T2); % Insert ad-hoc line broadening
                            else
                                fid{idxMolecule}(idxTime) = fid{idxMolecule}(idxTime) + trace(spins.molecule(idxMolecule).spin(idxSpin).rho*Ixy_tot)*exp(-curTime/T2); % Insert ad-hoc line broadening
                            end
                    end
                    % =============================================
                    % Compute gradient field & inhomogeneity offset
                    % =============================================
                    
                    %         __________________ ppm ___________________   __ T _  
                    %        /                                          \ /      \ 
                    offset = spins.molecule(idxMolecule).spin(idxSpin).B0*spins.B0; % in uT
                                   
                    %  pulse.Gx is in mT/m, Iz_tot_grd is in kHz/mT, and position vectors are in mm.
                    HGrad = 2*pi*(pulse.Gx(idxTime)*spins.molecule(idxMolecule).spin(idxSpin).r(1) ...
                                + pulse.Gy(idxTime)*spins.molecule(idxMolecule).spin(idxSpin).r(2) ...
                                + pulse.Gz(idxTime)*spins.molecule(idxMolecule).spin(idxSpin).r(3) ...
                                + offset)*0.001*Iz_tot_grd;  % In rad*kHz
                    % ==================
                    % Compute propagator
                    % ==================
                    if ((acqType~=0) && isSecular) 
                        P = (-1)^(isLHRot-1)*1i*(H0 + HGrad)*dt; % No RF. dt is in ms. HGrad, H0 is in rad*kHz.
                        U = diag(exp(diag(P)));
                    else
                        U = expm((-1)^(isLHRot-1)*1i*(H0 + HGrad + B1Scaling*H_RF{idxTime})*dt); 
                    end
                    % =========
                    % Propagate
                    % =========
                    if iscell(spins.molecule(idxMolecule).spin(idxSpin).rho)
                        spins.molecule(idxMolecule).spin(idxSpin).rho{idxSys} = U*spins.molecule(idxMolecule).spin(idxSpin).rho{idxSys}*U';
                    else
                        spins.molecule(idxMolecule).spin(idxSpin).rho = U*spins.molecule(idxMolecule).spin(idxSpin).rho*U';
                    end
                end % Loop over time
            end % Loop over spins
        end % Loop over subsystems
    end % Loop over molecules
end
