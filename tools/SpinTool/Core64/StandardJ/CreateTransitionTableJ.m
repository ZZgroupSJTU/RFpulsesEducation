function [TT, data, xAxis] = CreateTransitionTableJ(spins, varargin)
% CreateTransitionTableJ  Calculates spin system transitions' amplitudes and phases.
%   TT = CreateTransitionTableJ(spins)
%   Returns an Nx2 matrix, TT, such that TT(i,1) is the transition's
%   frequency (in Hz), TT(i,2) its (complex) amplitude (=amplitude & phase). 
%   It is assumed no pulse is applied during the creation of the table.
%
%   [TT,data,xAxis] = CreateTransitionTableJ(spins, 'paramName', paramValue, ...)
%   Parameter name-value pairs enabling additional functionality:
%     AcqType        A string which can be 'combined' or 'separate'. If set
%                    to combined, all molecules' transitions would be stored 
%                    in the same table, while separated would generate a cell 
%                    array of transition tables, one for each of the molecules 
%                    in the supplied spins. Default: 'combined'.
%     Domain         'freq' (default) or 'time'. Set the domain in which the 
%                    data (either FID or spectrum) is generated. 
%     T2             T2 decoherence time in milliseconds. If both T2 and
%                    FWHM are provided, FWHM takes precedence and T2 is 
%                    ignored. Default: empty set (i.e. use FWHM).
%                    Notes: 
%                    [1] If Gaussian lineshapes are employed, T2 takes on
%                    the meaning of the sigma parameter of the Gaussian
%                    distribution.
%                    [2] if Voigt lineshapes are employed, then T2
%                    should contain two elements, the first being the
%                    Lorentzian T2 and the second being the Gaussian
%                    sigma parameter.
%     isSeparateSys  Logical. Default: false. If set to true, each cell
%                    in the output data will contain a cell-array of the
%                    FIDs of each of the sub-systems of the molecule.
%                    For example, if the 3rd molecule is NAA, then
%                    data{3} will have two sub-cells, data{3}{1} and 
%                    data{3}{2}, corresponding to the singlet and 
%                    multiplet parts of NAA (the number of sub-systems will
%                    depend on how these were defined in the spin structure;
%                    consult the function SpinsJAddMolecule for examples).
%     FWHM           Full width at half max, in Hz. If omitted, it will be
%                    calculated based on T2 assuming a Lorentzian, via
%                    FWHM=1/(pi*T2). Default: 2 Hz.
%                    Note that Voigt lineshapes require TWO linewidths,
%                    so FWHM is a two-element vector (lorentzian, gaussian).
%     isCorrectDC    false by default. If set to true, any DC offset of the
%                    resulting wavefunctions in the spectral domain will
%                    be removed by subtracting the value of the farthest 
%                    point (with the highest frequency) in the frequency
%                    domain.
%     TimeAxis       Vector of time points at which to sample the FID, in ms. 
%                    Used if Domain='time'. If an empty set is used, then
%                    a 512 ms vector with 512 points is used. Default: 
%                    empty set. 
%     FreqAxis       Vector of frequency points (in kHz) at which to sample 
%                    the spectrum. Used if Domain='freq'. If an empty set 
%                    is used, FreqAxis is created with a spectral width of
%                    1000 points between +-0.5 kHz. Default: empty set. 
%     dV             Volume element per spin. This is used to normalize the
%                    signal if the number of spins changes but the 
%                    size of the sample doesn't (useful for visual 
%                    comparisons). Default: 1. 
%     LineShape      'lorentzian', 'gaussian', or 'voigt'. Default is
%                    'lorentzian'. Note that voigt lineshapes require a 
%                    TWO-element FWHM, since you must specify both the 
%                    linewidths of the original Lorentzian lineshape AND 
%                    the Gaussian convolution kernel.
%     isOutputCell   Logical (default: true). If true, data will be a
%                    cell array, one for each metabolite/subsystem.
%                    If false, data will be a matrix where each column
%                    is a separate metabolite (or just one where all
%                    metabolites are summed if acqType='combined').
%
%   data is either the FID or the spectrum of the molecules (depending on
%   the value of 'Domain'), and is in a cell array format. xAxis is the 
%   x-Axis (either in time or frequency domain) used to construct the data.

q = inputParser;
q.addParameter('AcqType', 'combined', @(x) ismember(lower(x), {'combined', 'separate'}));
q.addParameter('TimeAxis', [], @(x) isempty(x) || isvector(x));
q.addParameter('FreqAxis', [], @(x) isempty(x) || isvector(x));
q.addParameter('isCorrectDC', false, @(x) islogical(x));
q.addParameter('isSeparateSys', false, @(x) islogical(x));
q.addParameter('T2', [], @(x) isnumeric(x));
q.addParameter('FWHM', 2, @(x) isnumeric(x));
q.addParameter('Domain', 'freq', @(x) isempty(x) || ismember(lower(x), {'freq', 'time'}));
q.addParameter('LineShape', 'lorentzian', @(x) ismember(lower(x), {'lorentzian', 'gaussian', 'voigt'}));
q.addParameter('dV', 1, @(x) isnumeric(x));
q.addParameter('isOutputCell', true, @(x) islogical(x));
q.parse(varargin{:});
inputParams = q.Results;

isSecular = spins.isSecular;
numMolecules = numel(spins.molecule);
acqType = inputParams.AcqType;

% =========================================================================
% Initialize transition table
% =========================================================================
switch lower(acqType)
    case 'combined'
        TT = [];
    case 'separate'
            for idxMolecule=1:numMolecules
                if inputParams.isSeparateSys
                    if iscell(spins.molecule(idxMolecule).csVec)
                        numSubSystems = numel(spins.molecule(idxMolecule).csVec);
                    else
                        numSubSystems = 1;
                    end
                    for idxSys=1:numSubSystems
                        TT{idxMolecule}{idxSys} = [];
                    end
                else
                    TT{idxMolecule} = [];
                end
            end
end

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
    numSpins = numel(spins.molecule(idxMolecule).spin);
    for idxSys=1:numSubSystems
        numNuclei = numel(nucleus{idxSys});
        [gmRatio, spin] = GetGyromagneticRatio(nucleus{idxSys});
        numDims = prod(spin*2+1);
        % Step I: Prepare operators
        Ix_tot = zeros(numDims, numDims);
        Iy_tot = zeros(numDims, numDims);
        for k=1:numNuclei
            Ix_tot = Ix_tot + IxN(k,numNuclei,spin);
            Iy_tot = Iy_tot + IyN(k,numNuclei,spin);
        end
        Ixy_tot = Ix_tot + 1i*Iy_tot;

        % ========================================================================
        % Compute CS Hamiltonian (sans B0 offsets)
        % ========================================================================

        Hcs = zeros(numDims, numDims);
        for p=1:numNuclei
            cs = (csVec{idxSys}(p) - spins.csCenter)*spins.B0*gmRatio(p)/1000; % in kHz
            Hcs = Hcs - 2*pi*cs*IzN(p,numNuclei,spin);
        end
        
        % ========================================================================
        % Compute J-Coupling Hamiltonian, in 2*pi*kHz
        % ========================================================================

        HJ  = zeros(numDims, numDims);
        for p=1:numNuclei
            for k=p+1:numNuclei
                HJ = HJ + 2*pi*JMatrix{idxSys}(p,k)*0.001*( ...
                     (1-isSecular)*(IxN(p,numNuclei,spin)*IxN(k,numNuclei,spin) + IyN(p,numNuclei,spin)*IyN(k,numNuclei,spin)) ...
                     + IzN(p,numNuclei,spin)*IzN(k,numNuclei,spin));
            end
        end
        HSpins = Hcs + HJ;

        % ========================================================================
        % Calculate the transition table.
        %
        % Let rho be the density matrix, H the Hamiltonian, O the observable. We can
        % diagonlize the Hamiltonian, H=V*D*V', such that 
        %   U(t)=expm(-1i*H*t) = V*expm(-1i*D*t)*V' = V*UD*V'
        % and hence
        %   s(t) = trace(rho(t)*O) 
        %        = trace(U*rho*U'*O)
        %        = trace(V*UD*V'*rho*V*UD'*V'*O)
        %        = trace(UD*rhoV*UD'*OV)           rhoV = V'*rho*V, OV = V'*O*V
        %        = sum_m         ( (UD*rhoV*UD'*OV)_{mm} )
        %        = sum_{m,n,p,q} ( UD(m,n)*rhoV(n,p)*UD'(p,q)*OV(q,m) )
        %        = sum_{m,p}     ( UD(m,m)*rhoV(m,p)*UD'(p,p)*OV(p,m) )
        %        = sum_{m,p}     exp(-1i*(Em-Ep)*t)*OV(p,m)*rhoV(m,p) 
        %        = sum_{m,p}     exp(-1i*(Em-Ep)*t)*(OV.')(m,p)*rhoV(m,p) 
        % ========================================================================

        [V, D] = eig(HSpins); % Eigenvalues are in rad*kHz
        D = real(diag(D)); % Eigenvalues of a Hermitian hamiltonian are always real.
        OV = V'*Ixy_tot*V;
        rhoV = zeros(numDims, numDims);
        for idxSpin=1:numSpins
            if iscell(spins.molecule(idxMolecule).spin(idxSpin).rho)
                rhoV = rhoV + spins.molecule(idxMolecule).spin(idxSpin).rho{idxSys};
            else
                rhoV = rhoV + spins.molecule(idxMolecule).spin(idxSpin).rho;
            end
        end
        rhoV = V'*rhoV*V;
        for idx1=1:numDims
            for idx2=1:numDims
                amp = rhoV(idx1,idx2)*OV(idx2,idx1);
                if abs(amp)>1e-10
                    switch lower(acqType)
                        case 'combined'
                            TT(end+1, :) = [(D(idx2)-D(idx1))/2/pi  amp];
                        case 'separate'
                            if inputParams.isSeparateSys
                                TT{idxMolecule}{idxSys}(end+1, :) = [(D(idx2)-D(idx1))/2/pi  amp];
                            else
                                TT{idxMolecule}(end+1, :) = [(D(idx2)-D(idx1))/2/pi  amp];
                            end
                    end
                end
            end
        end
    end
end

if nargout==1
    return
end

dV = inputParams.dV;

counterBasisFunction = 0;
switch lower(inputParams.Domain)
    case 'freq'
        % =================================================================
        % Generate acquired data in frequency domain
        % =================================================================
        if isempty(inputParams.FreqAxis)
            SW = 1;
            numPts = 1000;
            dv = SW/numPts;
            freqAxis = [-SW/2:dv:SW/2-dv];
        else
            freqAxis = inputParams.FreqAxis;
        end
        xAxis = freqAxis;
        numAcqPoints = numel(freqAxis);
        switch lower(acqType)
            case 'separate'
                % =========================================================
                % Output each molecule (and maybe sub-system) separately
                % =========================================================
                for idxMolecule=1:numMolecules
                    if inputParams.isSeparateSys
                        numSubSystems = numel(TT{idxMolecule});
                        for idxSys=1:numSubSystems
                            if inputParams.isOutputCell
                                data{idxMolecule}{idxSys} = zeros(1, numAcqPoints);
                            else
                                counterBasisFunction = counterBasisFunction+1;
                                data(:,counterBasisFunction) = zeros(numAcqPoints,1);
                            end
                            numLines = size(TT{idxMolecule}{idxSys},1);
                            for idxLine=1:numLines
                                amp = TT{idxMolecule}{idxSys}(idxLine, 2)*dV; 
                                freq = -TT{idxMolecule}{idxSys}(idxLine,1);
                                curLineshape = GenerateLineshape('lineShape', inputParams.LineShape, ...
                                                  'amp', amp*dV, ...
                                                  'freqAxis', freqAxis, ...
                                                  'centerFreq', freq, ...
                                                  'FWHM', inputParams.FWHM, ...
                                                  'T2', inputParams.T2);
                                if inputParams.isOutputCell
                                    data{idxMolecule}{idxSys} = data{idxMolecule}{idxSys} + curLineshape;
                                else
                                    data(:,counterBasisFunction) = data(:,counterBasisFunction) + curLineshape.';
                                end
                            end
                        end
                    else
                        if inputParams.isOutputCell
                            data{idxMolecule} = zeros(1, numAcqPoints);
                        else
                            data(:, idxMolecule) = zeros(numAcqPoints, 1);
                        end
                        numLines = size(TT{idxMolecule},1);
                        for idxLine=1:numLines
                            amp = TT{idxMolecule}(idxLine, 2); 
                            freq = -TT{idxMolecule}(idxLine,1);
                            curLineshape = GenerateLineshape('lineShape', inputParams.LineShape, ...
                                              'amp', amp*dV, ...
                                              'freqAxis', freqAxis, ...
                                              'centerFreq', freq, ...
                                              'FWHM', inputParams.FWHM, ...
                                              'T2', inputParams.T2);
                            if inputParams.isOutputCell
                                data{idxMolecule} = data{idxMolecule} + curLineshape;
                            else
                                data(:, idxMolecule) = data(:, idxMolecule) + curLineshape.';
                            end
                        end
                    end
                end
            case 'combined'
                % =========================================================
                % Create one spectrum comprised of all molecules
                % =========================================================
                numLines = size(TT,1);
                data = zeros(1, numAcqPoints);
                for idxLine=1:numLines
                    amp = TT(idxLine, 2); 
                    freq = -TT(idxLine,1);
                    data = data + GenerateLineshape('lineShape', inputParams.LineShape, ...
                                                    'amp', amp*dV, ...
                                                    'freqAxis', freqAxis, ...
                                                    'centerFreq', freq, ...
                                                    'FWHM', inputParams.FWHM, ...
                                                    'T2', inputParams.T2);
                end
        end
    case 'time'
        % =================================================================
        % Generate acquired data in time domain
        % =================================================================
        if isempty(inputParams.TimeAxis)
            acqTime = 512;
            numAcqPoints = 512;
            dt = acqTime/numAcqPoints;
            timeAxis = [0:dt:(numAcqPoints-1)*dt];
        else
            timeAxis = inputParams.TimeAxis;
            numAcqPoints = numel(timeAxis);
        end
        xAxis = timeAxis;
        switch lower(acqType)
            case 'separate'
                % =========================================================
                % Output each molecule (and maybe sub-system) separately
                % =========================================================
                for idxMolecule=1:numMolecules
                    if inputParams.isSeparateSys
                        numSubSystems = numel(TT{idxMolecule});
                        for idxSys=1:numSubSystems
                            if inputParams.isOutputCell
                                data{idxMolecule}{idxSys} = zeros(1, numAcqPoints);
                            else
                                counterBasisFunction = counterBasisFunction+1;
                                data(:,counterBasisFunction) = zeros(numAcqPoints,1);
                            end
                            numLines = size(TT{idxMolecule}{idxSys},1);
                            for idxLine=1:numLines
                                freq = -TT{idxMolecule}{idxSys}(idxLine,1);
                                amp  = TT{idxMolecule}{idxSys}(idxLine,2); 
                                curLineshape = GenerateLineshape('lineShape', inputParams.LineShape, ...
                                                  'amp', amp*dV, ...
                                                  'timeAxis', timeAxis, ...
                                                  'centerFreq', freq, ...
                                                  'FWHM', inputParams.FWHM, ...
                                                  'T2', inputParams.T2);
                                if inputParams.isOutputCell
                                    data{idxMolecule}{idxSys} = data{idxMolecule}{idxSys} + curLineshape;
                                else
                                    data(:,counterBasisFunction) = data(:,counterBasisFunction) + curLineshape.';
                                end
                                                  
                            end
                        end
                    else
                        if inputParams.isOutputCell
                            data{idxMolecule} = zeros(1, numAcqPoints);
                        else
                            data(:, idxMolecule) = zeros(numAcqPoints, 1);
                        end
                        numLines = size(TT{idxMolecule},1);
                        for idxLine=1:numLines
                            freq = -TT{idxMolecule}(idxLine,1);
                            amp  = TT{idxMolecule}(idxLine, 2); 
                            curLineshape = GenerateLineshape('lineShape', inputParams.LineShape, ...
                                              'amp', amp*dV, ...
                                              'timeAxis', timeAxis, ...
                                              'centerFreq', freq, ...
                                              'FWHM', inputParams.FWHM, ...
                                              'T2', inputParams.T2);
                            if inputParams.isOutputCell
                                data{idxMolecule} = data{idxMolecule} + curLineshape;
                            else
                                data(:, idxMolecule) = data(:, idxMolecule) + curLineshape.';
                            end
                        end
                    end
                end
            case 'combined'
                % =========================================================
                % Create one spectrum comprised of all molecules
                % =========================================================
                numLines = size(TT,1);
                data = zeros(numAcqPoints, 1);
                for idxLine=1:numLines
                    amp = TT(idxLine,2);
                    freq = -TT(idxLine,1);
                    data = data + ...
                        GenerateLineshape('lineShape', inputParams.LineShape, ...
                                          'amp', amp*dV, ...
                                          'timeAxis', timeAxis, ...
                                          'centerFreq', freq, ...
                                          'FWHM', inputParams.FWHM, ...
                                          'T2', inputParams.T2).';
                end
        end
end

