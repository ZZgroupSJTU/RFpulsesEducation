function [spins,fidCellArray]=ApplySequence3D(spins, seq)
% [spins,fidCellArray]=ApplySequence3D(spins, seq)
%
% Applies a sequence seq to a 3D spin structure. Outputs the spin structure at
% the end of the sequence. Also returned is a cell array of FIDs acquired
% throughout the sequence.
%
% For the structure and possible entries of seq, consult ApplySequence.
% Note that this function operates on a simplified 3D spin structure, 
% created with InitSpins3D, and uses the appropriate ApplyPulse3D.
%
% IMPORTANT NOTE:
% The Spins3D has an additional command 'SumMxy', which takes no parameters, 
% and simply sums the transverse magnetization (without affecting it) when
% applied. 

numElements = numel(seq);
fidCellArray = {};
gmr = GetGyromagneticRatio('1h');

for idxEl=1:numElements
    curElement = seq{idxEl};
    if iscell(curElement)
        numVars = numel(curElement);
        cmd = curElement{1};
        M = squeeze(sum(sum(sum(spins.M,1),2),3)).';
        switch lower(cmd)
            case 'pulse'
                switch numVars
                    case 2 % {'pulse', pulse}
                        spins = ApplyPulse3D(spins, curElement{2});
                    case 3 % {'pulse', pulse, pulsePhase}
                        spins = ApplyPulse3DCycle(spins, curElement{2}, curElement{3}); 
                    case 4 % {'pulse', pulse, pulsePhase, addCoeffs}
                        spins = ApplyPulse3DCycle(spins, curElement{2}, curElement{3}, curElement{4}); 
                    case 5 % {'pulse', pulse, pulsePhase, addCoeffs, cycleType}
                        spins = ApplyPulse3DCycle(spins, curElement{2}, curElement{3}, curElement{4}, curElement{5}); 
                    case 6 % {'pulse', pulse, pulsePhase, addCoeffs, cycleType, B1}
                        spins = ApplyPulse3DCycle(spins, curElement{2}, curElement{3}, curElement{4}, curElement{5}, curElement{6}); 
                    otherwise
                        error('Number of variables %d in pulse element must be 2, 3 or 4.', numVars);
                end
            case 'delay'
                switch numVars
                    case 2 % {'delay', delay (ms)}
                        spins = Delay3D(spins, curElement{2});
                    case 3 % {'delay', delay (ms), numSteps}
                        pulseDelay = PulseCreateZero(curElement{2}, curElement{3});
                        spins = ApplyPulse3D(spins, pulseDelay);
                    otherwise
                        error('Number of variables %d in delay element must be 2 or 3.', numVars);
                end
            case 'hard'
                numSteps = 1;
                switch numVars
                    case 2 % {'hard', ang}
                        pulsePhase = 0; % deg.
                        pulseAngle = curElement{2}; % deg.
                        pulseDuration = 0.1; % us
                        spins = ApplyPulse3DHard(spins, pulseDuration, pulseAngle, pulsePhase);
                    case 3 % {'hard', ang, ph}
                        pulsePhase = curElement{3}; % deg.
                        pulseAngle = curElement{2}; % deg.
                        pulseDuration = 0.1; % us
                        pulseHard = PulseCreateHard(pulseDuration, pulseAngle, 0.0, numSteps);
                        spins = ApplyPulse3DCycle(spins, pulseHard, pulsePhase);  
                    case 4 % {'hard', ang, ph, addCoeffs}  
                        addCoeffs = curElement{4}; 
                        pulsePhase = curElement{3}; % deg.
                        pulseAngle = curElement{2}; % deg.
                        pulseDuration = 0.1; % us
                        pulseHard = PulseCreateHard(pulseDuration, pulseAngle, 0.0, numSteps);
                        spins = ApplyPulse3DCycle(spins, pulseHard, pulsePhase, addCoeffs);  
                    case 5 % {'hard', ang, ph, addCoeffs, cycleType}  
                        cycleType = curElement{5};
                        addCoeffs = curElement{4}; 
                        pulsePhase = curElement{3}; % deg.
                        pulseAngle = curElement{2}; % deg.
                        pulseDuration = 0.1; % us
                        pulseHard = PulseCreateHard(pulseDuration, pulseAngle, 0.0, numSteps);
                        spins = ApplyPulse3DCycle(spins, pulseHard, pulsePhase, addCoeffs, cycleType);  
                    case 6 % {'hard', ang, ph, addCoeffs, cycleType, B1}  
                        B1 = curElement{6};
                        cycleType = curElement{5};
                        addCoeffs = curElement{4}; 
                        pulsePhase = curElement{3}; % deg.
                        pulseAngle = curElement{2}; % deg.
                        pulseDuration = 0.1; % us
                        pulseHard = PulseCreateHard(pulseDuration, pulseAngle, 0.0, numSteps);
                        spins = ApplyPulse3DCycle(spins, pulseHard, pulsePhase, addCoeffs, cycleType, B1);  
                    otherwise
                        error('Number of variables %d in hard pulse sequence command should be 2,3,4.', numVars);
                end
            case 'rect'
                B1 = [];
                switch numVars
                    case 4 % {'rect', ang, ph, d} 
                        numSteps = 1;
                        numCycles = numel(curElement{3});
                        addCoeffs = ones(1,numCycles)/numCycles;
                        cycleType = 'normal';
                    case 5 % {'rect', ang, ph, d, numSteps}   
                        numSteps = curElement{5};
                        numCycles = numel(curElement{3});
                        addCoeffs = ones(1,numCycles)/numCycles;
                        cycleType = 'normal';
                    case 6 % {'rect', ang, ph, d, numSteps, addCoeffs}   
                        numSteps = curElement{5};
                        addCoeffs = curElement{6};
                        cycleType = 'normal';
                    case 7 % {'rect', ang, ph, d, numSteps, addCoeffs, cycleType}   
                        numSteps = curElement{5};
                        addCoeffs = curElement{6};
                        cycleType = curElement{7};
                    case 8 % {'rect', ang, ph, d, numSteps, addCoeffs, cycleType, B1}   
                        numSteps = curElement{5};
                        addCoeffs = curElement{6};
                        cycleType = curElement{7};
                        B1 = curElement{8};
                    otherwise
                        error('Number of variables %d in rect pulse should be 4, 5, 6, 7 or 8.', numVars);
                end
                pulseDuration = curElement{4}; % ms
                pulsePhase = curElement{3}; % deg.
                flipAngle = curElement{2}/360; % fraction of 360
                peakB1 = flipAngle/pulseDuration; % kHz
                rectPulse = PulseCreateConst(pulseDuration, numSteps, peakB1, 0);
                spins = ApplyPulse3DCycle(spins, rectPulse, pulsePhase, addCoeffs, cycleType, B1);
            case 'killmxy'
                spins = PurgePerfect3D(spins);
            case 'killmz'
                for idx=1:numel(spins)
                    spins(idx).M(3) = 0;
                end
            case 'thermal'
                for idx=1:numel(spins)
                    spins(idx).M = [0; 0; spins(idx).M0];
                end
            case 'purge'
                Gx = curElement{2}*gmr/1000; % Convert to kHz/mm
                Gy = curElement{3}*gmr/1000; % Convert to kHz/mm
                Gz = curElement{4}*gmr/1000; % Convert to kHz/mm
                tp = curElement{5};
                spins = Purge3D(spins, Gx, Gy, Gz, tp);
            case 'purgemoment'
                kx = curElement{2};
                ky = curElement{3};
                kz = curElement{4};
                tp = curElement{5};
                spins = PurgeMoment3D(spins, kx, ky, kz, tp);
            case 'summxy'
                summedM = squeeze(sum(sum(sum(spins.M, 1),2),3));
                fidCellArray{end+1} = summedM(1) + 1i*summedM(2);
            case 'acquire'
                if ispulse(curElement{2})  % {'acquire', pulse}
                    [spins, fid] = Acquire3D(spins, curElement{2});
                    fidCellArray = [fidCellArray, fid];
                else % {'acquire', Nt, SW, Gx, Gy, Gz}
                    Gx = curElement{4}*gmr/1000; % Convert to kHz/mm
                    Gy = curElement{5}*gmr/1000; % Convert to kHz/mm
                    Gz = curElement{6}*gmr/1000; % Convert to kHz/mm
                    numAcqPoints = curElement{2};
                    SW = curElement{3};
                    dt = 1/SW;
                    acqTime = dt*numAcqPoints;
                    [spins, fid] = Acquire3DGrad(spins, acqTime, numAcqPoints, Gx, Gy, Gz);
                    fidCellArray = [fidCellArray, fid];
                end
            otherwise
                fprintf('Unrecognized command in sequence: %s \n', cmd);
        end
    else % A pulse
        spins = ApplyPulse3D(spins, curElement);
    end
end