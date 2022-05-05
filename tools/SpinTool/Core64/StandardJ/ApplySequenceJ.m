function [spins,fidCellArray]=ApplySequenceJ(spins, seq, isFeedbackOn)
% [spins,fidCellArray]=ApplySequenceJ(spins, seq, isFeedbackOn)
%
% Applies a sequence seq to a J-coupled spin structure. Outputs the spin 
% structure at the end of the sequence. Also returned is a cell array of 
% FIDs acquired throughout the sequence.
%
% isFeedbackOn is a logical (boolean) value, false by default. If true,
% feedback is provided for each step of the simulation.
%
% A Sequence is a cell array containing structures/commands that are acted
% out serially on the given spin structure. The possible cell elements are
% either a pulse structure or a string command. If a pulse structure is
% encountered it will be applied to the spins. The following string
% commands are possible and are case insensitive:
%
%   {'pulse', p, [affectedNuclei], [Phase Cycle], [Addition Coefficients], [freqRange]}
%      Applies a pulse p to the spins. Optional arguments are:
%        affectedNuclei: A 1xM cell array, where M = number of molecules.
%                        Each element in the cell array is a 1xN logical 
%                        array, with N = number of nuclei in that molecule.
%                        Each element in the array signifies whether the
%                        pulse will be applied to that particular nucleus.
%                        Omit, or set to [], to affect all nuclei.
%        Phase Cycle:    An array of pulse phases, in degrees. The 
%                        pulse will be applied with varying phases. Omit,
%                        or set to [], to use a single cycle with 0 phase.
%        Add. Coeff.:    Array with the same number of elements as Phase
%                        Cycle, indicating the weights of each path
%                        (complex weights possible!). If omitted, 
%                        each path will be given equal weight, with the
%                        sum of all weights = 1.
%        freqRange       A 1x2 vector of min ppm, max ppm which the pulse
%                        will affect. This is a crude way of simulating
%                        frequency selective pulses, and is often more 
%                        intuitive than setting affectedNuclei directly.
%
%   {'delay', d}           
%      Applies a delay for d milliseconds.
%
%   {'hard', ang, ph, [affectedNuclei], [freqRange]}
%      Applies a hard pulse (0.1 microsecs) with a flip angle ang (deg) 
%      and phase ph (deg). 
%       
%   {'rect', ang, ph, d, [affectedNuclei]}   
%      Applies a rect pulse with a flip angle ang (deg) phase ph (deg) and 
%      duration d (ms). 
%      For a description of affectedSpins, see the entry for 'pulse'.
% 
%   {'const', ang, ph, [numSteps], [affectedNuclei]}
%      Applies a constant pulse designed to achieve a certain flip angle
%      and (deg) with a certain pulse phase (deg), and a certain number of 
%      steps (numSteps). 
%      For a description of affectedNuclei, see the entry for 'pulse'.
% 
%   {'gaussian', ang, ph, centerFreq, FWHM, SW, [affectedNuclei]}
%      Applies a frequency selected Gaussian pulse to the spectrum with
%      a given FWHM and SW in Hz.
%
%   {'purge', Gx, Gy, Gz, d}
%      Applies a delay d (ms) with gradients Gx, Gy, Gz along the x, y and
%      z axes (in mT/m)
%
%   {'purgemoment', kx, ky, kz, d, nucleus}
%      Applies a delay d (ms) with gradient moments kx, ky, kz along the
%      x, y, z axes (in m^(-1)). The gradients are calibrated for the given
%      nucleus: '1h', '13c', '31p', etc ... depending on the gyromagnetic
%      ratio.
%  
%   {'killmxy', [affectedNuclei]}
%      Kills all magnetization in the xy plane.
%      For a description of affectedSpins, see the entry for 'pulse'.
%
%   {'acquire', Nt, SW, Gx, Gy, Gz}
%      Acquires an FID with Nt points, SW spectral width (kHz) and 
%      gradients Gx, Gy and Gz along the x, y, z axes (in mT/m).
%
%   Example:
%   p = PulseCreateHSn(1, 5, 5.3, 7, 3, 512); % Create hyperbolic secant
%   seq = {{'pulse', p}, ...
%          {'delay', 10}, ...
%          {'acquire', 1000, 1, 0, 0, 0}}
%   spins = InitSpinsRelax(0, 500, 100, [0; 0; 1], 1e6, 1e6, 1, 0, 1);
%   [spins, fidCellArray] = ApplySequenceJ(spins, seq) % fidCellArray will
%                                                     % only have 1 element

if nargin<3, isFeedbackOn = false; end

numElements = numel(seq);
fidCellArray = {};

if isFeedbackOn, fprintf('Simulation started.\n'); end

for idxEl=1:numElements
    curElement = seq{idxEl};
    if iscell(curElement)
        cmd = curElement{1};
        switch lower(cmd)
            case 'delay'
                spins = DelayJ(spins, curElement{2});
                if isFeedbackOn, fprintf('  Delay complete.\n'); end
            case 'hard' % {'hard', ang, ph, [affectedNuclei], [freqRange], [addition coefficients]}
                numParams = numel(curElement); 
                if numParams>6, error('Problem executing hard pulse in ApplySequenceJ: %d parameters supplied, but hard pulses can only accept a maximum of 6', numParams); end
                tiltAngle = curElement{2}; % Deg.
                pulsePhase = curElement{3}; % Deg.
                pulseDuration = 0.0001; % ms
                if numParams<4, affectedNuclei = []; else, affectedNuclei = curElement{4}; end
                if numParams<5, freqRange = []; else, freqRange = curElement{5}; end
                numCycles = numel(pulsePhase);
                if numParams<6, addCoeff = ones(1, numCycles)/numCycles; end
                pulse = PulseCreateHard(pulseDuration, tiltAngle, 0.0);
                spins = ApplyPulseCycleJ(spins, pulse, affectedNuclei, pulsePhase, addCoeff, freqRange);
                if isFeedbackOn, fprintf('  Hard pulse complete.\n'); end
            case 'rect' % {'rect', ang, ph, d, [affectedNuclei], [freqRange]}   
                numParams = numel(curElement); 
                tiltAngle = curElement{2}; % Deg.
                pulsePhase = curElement{3}; % Deg.
                pulseDuration = curElement{4}; % ms
                if numParams<5, affectedNuclei = []; else, affectedNuclei = curElement{5}; end
                if numParams<6, freqRange = []; else, freqRange = curElement{6}; end
                spins = ApplyHardPulseJ(spins, pulseDuration, tiltAngle, pulsePhase, affectedNuclei, freqRange);
                if isFeedbackOn, fprintf('  Rect pulse complete.\n'); end
            case 'gaussian'  %  {'gaussian', ang, ph, centerFreq, FWHM, SW, [affectedNuclei]}
                numParams = numel(curElement); 
                tiltAngle = curElement{2}; % Deg.
                phaseCycle = curElement{3}; % Deg.
                centerFreq = curElement{4}; % Hz
                FWHM = curElement{5}; % Hz
                SW = curElement{6}; % Hz
                if numParams<7, affectedNuclei = []; else affectedNuclei = curElement{7}; end
                pulse = PulseCreateGaussian(centerFreq*0.001, FWHM*0.001, SW*0.001, tiltAngle);
                spins = ApplyPulseCycleJ(spins, pulse, affectedNuclei, phaseCycle, addCoeffs);
                if isFeedbackOn, fprintf('  Gaussian pulse complete.\n'); end
            case 'killmxy' % {'killmxy', [affectedNuclei]}
                error('killmxy Not yet implemented.');
                numParams = numel(curElement); 
                if numParams<2
                    spins = PurgePerfectJ(spins,0.0000001);
                else
                    affectedNuclei = curElement{4};
                    spins = PurgePerfectJ(spins,0.0000001,affectedNuclei);
                end
            case 'purge'
                Gx = curElement{2}; % mT/m
                Gy = curElement{3}; % mT/m
                Gz = curElement{4}; % mT/m
                tp = curElement{5};
                spins = DelayJ(spins, tp, Gx, Gy, Gz); 
                if isFeedbackOn, fprintf('  Purge complete.\n'); end
            case 'purgemoment' % {'purgemoment', kx, ky, kz, d, [nucleus]}
                numParams = numel(curElement); 
                kx = curElement{2};
                ky = curElement{3};
                kz = curElement{4};
                tp = curElement{5};
                if numParams<6, nucleus='1h'; else, nucleus = curElement{6}; end
                Gx = 0.001*kx/tp*GetGyromagneticRatio('1h')/GetGyromagneticRatio(nucleus); % kHz/mm
                Gy = 0.001*ky/tp*GetGyromagneticRatio('1h')/GetGyromagneticRatio(nucleus); % kHz/mm
                Gz = 0.001*kz/tp*GetGyromagneticRatio('1h')/GetGyromagneticRatio(nucleus); % kHz/mm
                % Gx = kx/(tp*GetGyromagneticRatio(nucleus)); 
                % Gy = ky/(tp*GetGyromagneticRatio(nucleus));
                % Gz = kz/(tp*GetGyromagneticRatio(nucleus));
                spins = DelayJ(spins, tp, Gx, Gy, Gz);
                if isFeedbackOn, fprintf('  PurgeMoment complete.\n'); end
            case 'acquire' % {'acquire', Nt, SW, Gx, Gy, Gz}
                Gx = curElement{4}; % mT/m
                Gy = curElement{5}; % mT/m
                Gz = curElement{6}; % mT/m
                numAcqPoints = curElement{2};
                SW = curElement{3}; % kHz
                acqTime = 1/SW*numAcqPoints;
                [fid, spins] = AcquireJ(spins, acqTime, numAcqPoints, Gx, Gy, Gz);
                fidCellArray = [fidCellArray, fid];
                if isFeedbackOn, fprintf('  Acquire complete.\n'); end
            case 'pulse' % {'pulse', p, [affectedNuclei], [Phase cycle], [Addition coefficients], [freqRange]}
                numParams = numel(curElement); 
                pulse = curElement{2};
                if numParams<3, affectedNuclei = []; else, affectedNuclei = curElement{3}; end
                if numParams<4, phaseCycle = 0; else, phaseCycle = curElement{4}; end
                if numParams<5, addCoeffs = []; else, addCoeffs = curElement{5}; end
                if numParams<6, freqRange = []; else, freqRange = curElement{6}; end
                spins = ApplyPulseCycleJ(spins, pulse, affectedNuclei, phaseCycle, addCoeffs, freqRange);
                if isFeedbackOn, fprintf('  Pulse complete.\n'); end
            otherwise
                fprintf('Unrecognized command in sequence: %s \n', cmd);
        end
    else % Unknown - kill it with fire!!!
        error('Sequence entry unrecognized.');
    end
end