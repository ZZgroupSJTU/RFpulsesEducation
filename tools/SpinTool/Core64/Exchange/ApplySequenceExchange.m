function [spins,fidCellArray]=ApplySequenceExchange(spins, seq)
% [spins,fidCellArray]=ApplySequenceExchange(spins, seq)
%
% Applies a sequence seq to a spin structure. Outputs the spin structure at
% the end of the sequence. Also returned is a cell array of FIDs acquired
% throughout the sequence.
%
% The spins are in a structure of the form (N: number of spins per position):
%   spins(i).r      [x; y; z] in mm
%   spins(i).M      3N*1 vector of the form [Mx1; My1; Mz1; Mx2; My2;
%                   Mz2; ... ]
%   spins(i).cs     N*1 vector of chemical shifts of spins, in kHz
%   spins(i).T1     N*1 vector of T1s of spins, ms
%   spins(i).T2     N*1 vector of T2s of spins, ms
%   spins(i).M0     N*1 vector of equilibrium magnetization values
%   spins(i).B1     Scales RF at the position
%   spins(i).RS     Receiver sensitivity scaling at position r
%   spins(i).K      N*N exchange matrix, with exchange constants in Hz
%
% A Sequence is a cell array containing structures/commands that are acted
% out serially on the given spin structure. The possible cell elements are
% either a pulse structure or a string command. If a pulse structure is
% encountered it will be applied to the spins. The following string
% commands are possible and are case insensitive:
%
%   {'delay', d, [numSteps]}
%   Applies a delay for d milliseconds and numSteps number of steps 
%   (numSteps defaults to 1).
%
%   {'hard', ang, ph}
%   Applies a hard pulse (0.1 microsecs) with a flip angle ang (deg) 
%   and phase ph (deg)
%
%   {'rect', ang, ph, d, [centerFreq], [numSteps]}    
%   Applies a rect pulse with a flip angle ang (deg) phase ph (deg) and 
%   duration d (ms). One can also set the center frequency (kHz) and
%   specify the number of steps (defaults to 1).
% 
%   {'rect', ang, ph, d, numSteps}   
%   Same as 'rect' but here a number of steps can be specified. 
% 
%   {'purge', Gx, Gy, Gz, d}
%   Applies a delay d (ms) with gradients Gx, Gy, Gz along the x, y and
%   z axes (in mT/m)
%
%   {'purgemoment', kx, ky, kz, d}
%   Applies a delay d (ms) with gradient moments kx, ky, kz along the
%   x, y, z axes (in m^(-1))
%  
%   {'killmz'}
%   Kills all longitudinal magnetization. Not physical, but useful for
%   debugging.
%
%   {'pulse', pulse, [centerFreq]}
%   Applies a pulse object, having the following structure:
%     pulse.tp        In milliseconds
%     pulse.RFamp     1xN vector in kHz
%     pulse.RFphase   1xN vector in radians
%     pulse.Gx        1xN vector in kHz/mm
%     pulse.Gy        1xN vector in kHz/mm 
%     pulse.Gz        1xN vector in kHz/mm
%
%   {'thermal'}
%   Returns spins to their thermal equilibrium (M0) values, and kills
%   all transverse magnetization.
%
%   {'killmxy'}
%   Kills all magnetization in the xy plane
%
%   {'killmz'}
%   Kills all magnetization along the z-axis.
%
%   {'acquire', Nt, SW, Gx, Gy, Gz}
%   Acquires an FID with Nt points, SW spectral width (kHz) and gradients
%   Gx, Gy and Gz along the x, y, z axes (in mT/m)
%
%   {'acquire', pulse}
%   Acquires an FID while playing out time-dependent RF and gradient
%   shapes. 

numElements = numel(seq);
fidCellArray = {};
gmr = GetGyromagneticRatio('1h');

for idxEl=1:numElements
    curElement = seq{idxEl};
    numElements = numel(curElement);
    if iscell(curElement)
        cmd = curElement{1};
        switch lower(cmd)
            case 'delay' % {'delay', d, [numSteps]}
                pulseDelay = curElement{2};
                if numElements<3, numSteps = 1; else numSteps = curElement{3}; end
                pulseDelay = PulseCreateZero(pulseDelay, numSteps);
                spins = ApplyPulseExchange(spins, pulseDelay);
            case 'hard' % {'hard', ang, ph}      
                pulseHard = PulseCreateHard(0.1, curElement{2}, curElement{3}, 1);
                spins = ApplyPulseExchange(spins, pulseHard);
            case 'rect' % {'rect', ang, ph, d, [centerFreq], [numSteps]}    
                if numel(curElement)<5, centerFreq = 0; else centerFreq = curElement{5}; end
                if numel(curElement)>=6, numSteps = curElement{6}; else numSteps = 1; end
                pulseDuration = curElement{4}; % ms
                pulsePhase = curElement{3}/180*pi; % radians
                flipAngle = curElement{2}/360; % fraction of 360
                peakB1 = flipAngle/pulseDuration; % kHz
                pulseRect = PulseCreateConst(pulseDuration, numSteps, peakB1, pulsePhase);
                spins = ApplyPulseExchange(spins, pulseRect, centerFreq);
            case 'killmxy'
                for idx=1:numel(spins)
                    spins(idx).M(1:3:end) = 0;
                    spins(idx).M(2:3:end) = 0;
                end
            case 'killmz'
                for idx=1:numel(spins)
                    spins(idx).M(3:3:end) = 0;
                end
            case 'thermal'
                for idx=1:numel(spins)
                    spins(idx).M(1:3:end) = 0;
                    spins(idx).M(2:3:end) = 0;
                    if numel(spins(idx).M0)==1
                        numSites = round(spins(idx).M/3);
                        M0 = ones(numSites, 1)*M0;
                    end
                    spins(idx).M(3:3:end) = M0;
                end
            case 'purge' % {'purge', Gx, Gy, Gz, d}
                Gx = curElement{2}; % kHz/mm
                Gy = curElement{3}; % kHz/mm
                Gz = curElement{4}; % kHz/mm
                tp = curElement{5}; % ms
                pulse.tp      = tp;
                pulse.RFamp   = [0]; 
                pulse.RFphase = [0];
                pulse.Gx      = [Gx];
                pulse.Gy      = [Gy];
                pulse.Gz      = [Gz];
                spins = ApplyPulseExchange(spins,pulse);
            case 'purgemoment' % {'purgemoment', kx, ky, kz, d}
                kx = curElement{2}; % m^(-1)
                ky = curElement{3}; % m^(-1)
                kz = curElement{4}; % m^(-1)
                tp = curElement{5}; % ms
                Gx = 0.001*kx/tp; % kHz/mm
                Gy = 0.001*ky/tp; % kHz/mm
                Gz = 0.001*kz/tp; % kHz/mm
                pulse.tp      = tp;
                pulse.RFamp   = [0]; 
                pulse.RFphase = [0];
                pulse.Gx      = [Gx];
                pulse.Gy      = [Gy];
                pulse.Gz      = [Gz];
                spins = ApplyPulseExchange(spins,pulse);
            case 'pulse' % {'pulse', pulse [, centerFreq]}
                if ~ispulse(curElement{2})
                    error('Sequence ''pulse'' command encountered, but 2nd argument is not a pulse!'); 
                end
                if numel(curElement)<3
                    pulseOffset = 0;
                else
                    pulseOffset = curElement{3};
                end
                spins = ApplyPulseExchange(spins, curElement{2}, pulseOffset);
            case 'acquire'  
                if ispulse(curElement{2}) % {'acquire', pulse}
                    [spins, fid] = ApplyPulseExchange(spins, curElement{2});
                    fidCellArray = [fidCellArray, fid];
                else % {'acquire', Nt, SW, Gx, Gy, Gz} 
                    Gx            = curElement{4}; % kHz/mm
                    Gy            = curElement{5}; % kHz/mm
                    Gz            = curElement{6}; % kHz/mm
                    numAcqPoints  = curElement{2};
                    SW            = curElement{3}; % kHz
                    dt            = 1/SW; % ms
                    acqTime       = dt*numAcqPoints;
                    pulse.tp      = acqTime;
                    pulse.RFamp   = zeros(1,numAcqPoints);
                    pulse.RFphase = zeros(1,numAcqPoints);
                    pulse.Gx      = Gx.*ones(1,numAcqPoints);
                    pulse.Gy      = Gy.*ones(1,numAcqPoints);
                    pulse.Gz      = Gz.*ones(1,numAcqPoints);
                    [spins, fid]  = ApplyPulseExchange(spins, pulse);
                    fidCellArray  = [fidCellArray, fid];
                end
            otherwise
                fprintf('Unrecognized command in sequence: %s \n', cmd);
        end
    else % A pulse
        spins = ApplyPulseExchange(spins, curElement);
    end
end