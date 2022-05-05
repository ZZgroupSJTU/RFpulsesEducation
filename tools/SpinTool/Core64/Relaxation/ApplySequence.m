function [spins,fidCellArray]=ApplySequence(spins, seq)
% [spins,fidCellArray]=ApplySequence(spins, seq)
%
% Applies a sequence seq to a spin structure. Outputs the spin structure at
% the end of the sequence. Also returned is a cell array of FIDs acquired
% throughout the sequence.
%
% A Sequence is a cell array containing structures/commands that are acted
% out serially on the given spin structure. The possible cell elements are
% either a pulse structure or a string command. If a pulse structure is
% encountered it will be applied to the spins. The following string
% commands are possible and are case insensitive:
%
%   {'delay', d, [numSteps]}           
%   Applies a delay for d milliseconds. numSteps is an optional integer
%   parameter signifying the number of time steps in the delay.
%
%   {'pulse', pulse, [pulsePhase], [addCoeff], [cycleType], [B1]}
%   Applies an RF pulse given by the pulse object pulse. This can also
%   be replaced with just the pulse object. For information about
%   pulse objects, look at ApplyPulseRelax.c
%   pulsePhase is an optional 1xN vector of phases (deg.), such that
%   a phase cycle is carried out. 
%   addCoeff is an optional 1xN vector of addition coefficients for the
%   phase cycles (if omitted, phase cycles are just added up and normalized
%   by the number of cycles)
%   cycleType is a string which allows one to execute non-physical phase
%   cycles. Possible values (case insensitive):
%     'normal'      Phase cycle transformation R(phi) = Rz(phi)*R(0)*Rz(-phi)
%                   Here for example R(0)-R(pi/2)+R(pi)-R(3pi/2) yields the
%                   -1-->+1 and 1-->(-1) transitions ONLY.
%     'inverted'    Phase cycle transformation R(phi) = Rz(phi)*R(0)*Rz(+phi)
%                   Here for example R(0)+R(pi/2)+R(pi)+R(3pi/2) yields the
%                   -1-->+1, 1-->(-1) and 0-->0 transitions, so retains
%                   the longitudinal magnetization.
%   Some examples for common phase cycling schemes:
%     pulsePhase        addCoeff       cycleType    Purpose
%     0, 90, 180, 270   [1 -1 1 -1]/4  normal       Leave only 1-->(-1) and
%                                                   (-1)-->(+1) pathways.
%                                                   Kills off longitudinal
%                                                   magnetization!
%     0, 180            [1 1]/2        normal       Kills off excitation and
%                                                   storage pathways
%                                                   (0-->(+-1), (+-1)-->0)
%                                                   Leave z-magnetiation
%                                                   intact. Leaves xy
%                                                   magnetization intact.
%     0, 90, 180, 270  [1 1 1 1]/4     inverted     A "perfect 180": leaves
%                                                   only (+1)<-->(-1) and
%                                                   0-->0 transitions.
%   B1 is a possible fixed value of B1 for the pulse. This can be used, for
%   example, to approximate adiabatic pulses by using hard pulses and
%   fixing B1 inhomogeneity at 1.0. Omit or set to [] to use the true 
%   B1 of the spins.
%
%   {'hard', ang, [ph], [addCoeff], [cycleType], [B1]}
%   Applies a hard pulse (0.1 microsecs) with a flip angle ang (deg).
%   ph (deg.) is an optional 1xN vector of pulse phases for a phase rotation
%   scheme, and addCoeffs are optional addition coefficients. For more
%   information, see 'pulse'
%
%   {'rect', ang, ph, d, [numSteps], [addCoeffs], [cycleType], [B1]}   
%   Applies a rect pulse with a flip angle ang (deg) phase ph (deg) and 
%   duration d (ms). The optional parameter numSteps allows one to
%   specify the number of steps used for the pulse (1 step is the default).
%   If ph is a 1xN vector, a phase cycle will be carried out, and added
%   up. For phase cycling information, as well as the meaning of the 
%   addCoeffs 1xN vector, see the 'pulse' command.
% 
%   {'purge', Gx, Gy, Gz, d}
%   Applies a delay d (ms) with gradients Gx, Gy, Gz along the x, y and
%   z axes (in mT/m)
%
%   {'purgemoment', kx, ky, kz, d}
%   Applies a delay d (ms) with gradient moments kx, ky, kz along the
%   x, y, z axes (in m^(-1))
%
%   {'purgecsi', fov, step, numvoxels, d}
%   Applies a delay d (ms) with gradient moments which encode a particular
%   CSI step. The fov is a 1x3 vector in mm, and step is the number of the
%   current phase encoding step. For example: the 4th step of the 16-step
%   CSI encoding (16-voxels) with a 100 mm FOV along the z-axis is:
%   {'purgecsi', [0 0 100], [0 0 4], [0 0 16], 0.1}
%  
%   {'killmz'}
%   Kills all longitudinal magnetization. Not physical, but useful for
%   debugging.
%
%   {'thermal'}
%   Returns spins to their thermal equilibrium (M0) values, and kills
%   all transverse magnetization.
%
%   {'killmxy'}
%   Kills all magnetization in the xy plane
%
%   {'acquire', Nt, SW, Gx, Gy, Gz}
%   Acquires an FID with Nt points, SW spectral width (kHz) and gradients
%   Gx, Gy and Gz along the x, y, z axes (in mT/m)
%
%   {'acquireFOV', res, FOV, acqTime, axis}
%   Executes a readout element along a given axis. res is an integer.
%   FOV is in mm. acqTime is in ms. axis is a character: 'x', 'y', 'z'.
%   Note that FOV can be negative (to enforce negative gradient polarity).
%
%   {'acquireBW', res, FOV, BWPerPixel, axis}
%   Much like acquireFOV, only here the acquisition time is replaced with
%   the bandwidth per pixel (in Hz). Acquisition time is determined
%   indirectly from the BW per pixel and other parameters. BWPerPixel is
%   always positive!
%
%   {'acquire', pulse}
%   Acquires an FID while playing out time-dependent RF and gradient
%   shapes
%
%   Example:
%   p = PulseCreateHSn(1, 5, 5.3, 7, 3, 512); % Create hyperbolic secant
%   seq = {p, ...
%          {'delay', 10}, ...
%          {'acquire', 1000, 1, 0, 0, 0}}
%   spins = InitSpinsRelax(0, 500, 100, [0; 0; 1], 1e6, 1e6, 1, 0, 1);
%   [spins, fidCellArray] = ApplySequence(spins, seq) % fidCellArray will
%                                                     % only have 1 element

numElements = numel(seq);
fidCellArray = {};
gmr = GetGyromagneticRatio('1h');

for idxEl=1:numElements
    curElement = seq{idxEl};
    if iscell(curElement)
        numVars = numel(curElement);
        cmd = curElement{1};
        switch lower(cmd)
            case 'pulse'
                switch numVars
                    case 2 % {'pulse', pulse}
                        spins = ApplyPulseRelax(spins, curElement{2});
                    case 3 % {'pulse', pulse, pulsePhase}
                        spins = ApplyPulseCycle(spins, curElement{2}, curElement{3}); 
                    case 4 % {'pulse', pulse, pulsePhase, addCoeffs}
                        spins = ApplyPulseCycle(spins, curElement{2}, curElement{3}, curElement{4}); 
                    case 5 % {'pulse', pulse, pulsePhase, addCoeffs, cycleType}
                        spins = ApplyPulseCycle(spins, curElement{2}, curElement{3}, curElement{4}, curElement{5}); 
                    case 6 % {'pulse', pulse, pulsePhase, addCoeffs, cycleType, B1}
                        spins = ApplyPulseCycle(spins, curElement{2}, curElement{3}, curElement{4}, curElement{5}, curElement{6}); 
                    otherwise
                        error('Number of variables %d in pulse element must be between two and five.', numVars);
                end
            case 'delay'
                switch numVars
                    case 2 % {'delay', delay (ms)}
                        spins = DelayRelax(spins, curElement{2});
                    case 3 % {'delay', delay (ms), numSteps}
                        pulseDelay = PulseCreateZero(curElement{2}, curElement{3});
                        spins = ApplyPulseRelax(spins, pulseDelay);
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
                        spins = ApplyPulseHard(spins, pulseDuration, pulseAngle, pulsePhase);
                    case 3 % {'hard', ang, ph}
                        pulsePhase = curElement{3}; % deg.
                        pulseAngle = curElement{2}; % deg.
                        pulseDuration = 0.1; % us
                        pulseHard = PulseCreateHard(pulseDuration, pulseAngle, 0.0, numSteps);
                        spins = ApplyPulseCycle(spins, pulseHard, pulsePhase);  
                    case 4 % {'hard', ang, ph, addCoeffs}  
                        addCoeffs = curElement{4}; 
                        pulsePhase = curElement{3}; % deg.
                        pulseAngle = curElement{2}; % deg.
                        pulseDuration = 0.1; % us
                        pulseHard = PulseCreateHard(pulseDuration, pulseAngle, 0.0, numSteps);
                        spins = ApplyPulseCycle(spins, pulseHard, pulsePhase, addCoeffs);  
                    case 5 % {'hard', ang, ph, addCoeffs, cycleType}  
                        cycleType = curElement{5};
                        addCoeffs = curElement{4}; 
                        pulsePhase = curElement{3}; % deg.
                        pulseAngle = curElement{2}; % deg.
                        pulseDuration = 0.1; % us
                        pulseHard = PulseCreateHard(pulseDuration, pulseAngle, 0.0, numSteps);
                        spins = ApplyPulseCycle(spins, pulseHard, pulsePhase, addCoeffs, cycleType);  
                    case 6 % {'hard', ang, ph, addCoeffs, cycleType, B1}  
                        B1 = curElement{6};
                        cycleType = curElement{5};
                        addCoeffs = curElement{4}; 
                        pulsePhase = curElement{3}; % deg.
                        pulseAngle = curElement{2}; % deg.
                        pulseDuration = 0.1; % us
                        pulseHard = PulseCreateHard(pulseDuration, pulseAngle, 0.0, numSteps);
                        spins = ApplyPulseCycle(spins, pulseHard, pulsePhase, addCoeffs, cycleType, B1);  
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
                spins = ApplyPulseCycle(spins, rectPulse, pulsePhase, addCoeffs, cycleType, B1);
            case 'killmxy'
                spins = PurgePerfect(spins,0.0000001);
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
                spins = PurgeRelax(spins, Gx, Gy, Gz, tp);
            case 'purgemoment'
                kx = curElement{2};
                ky = curElement{3};
                kz = curElement{4};
                tp = curElement{5};
                spins = PurgeMoment(spins, kx, ky, kz, tp);
            case 'purgecsi'
                FOV = curElement{2}; % mm
                numVoxels = curElement{3};
                CSIStep = curElement{4}; 
                d = curElement{5}; % ms
                k = [0 0 0];
                for idxAxis=1:3
                    if FOV(idxAxis)>0 && numVoxels(idxAxis)>0
                        voxelSize = FOV(idxAxis)/numVoxels(idxAxis);
                        dk = 1/(FOV(idxAxis)*0.001); % m^(-1)
                        kMax = 1/(voxelSize*0.001); % m^(-1)
                        k(idxAxis) = -kMax/2 + CSIStep(idxAxis)*dk;
                    end
                end
                spins = PurgeMoment(spins, k(1), k(2), k(3), d);
            case 'acquire'
                if ispulse(curElement{2})  % {'acquire', pulse}
                    [spins, fid] = AcquirePulseRelax(spins, curElement{2});
                    fidCellArray = [fidCellArray, fid];
                else % {'acquire', Nt, SW, Gx, Gy, Gz}
                    Gx = curElement{4}*gmr/1000; % Convert to kHz/mm
                    Gy = curElement{5}*gmr/1000; % Convert to kHz/mm
                    Gz = curElement{6}*gmr/1000; % Convert to kHz/mm
                    numAcqPoints = curElement{2};
                    SW = curElement{3};
                    dt = 1/SW;
                    acqTime = dt*numAcqPoints;
                    [spins, fid] = AcquireGradRelax(spins, acqTime, numAcqPoints, Gx, Gy, Gz);
                    fidCellArray = [fidCellArray, fid];
                end
            case 'acquirebw' % {'acquireBW', res, FOV, BWPerPixel, axis}
                FOV = curElement{3}; % mm
                res = curElement{2}; % Integer
                BWPerPixel = curElement{4}*0.001; % kHz
                dx = FOV/res; % mm
                kMax = 1/dx; % 1/mm
                G = 1/(BWPerPixel*dx); % kHz/mm
                switch lower(curElement{5})
                    case 'x'
                        Gx = G; % kHz/mm 
                        Gy = 0; 
                        Gz = 0;
                    case 'y'
                        Gy = G; % kHz/mm
                        Gz = 0; 
                        Gx = 0;
                    case 'z'
                        Gz = G; % kHz/mm
                        Gx = 0; 
                        Gy = 0;
                    otherwise
                        error('Unrecognized axis %s in acquireFOV command', curElement{5});
                end
                acqTime = kMax/G; % ms (=kMax/Gradient)
                [spins, fid] = AcquireGradRelax(spins, acqTime, res, Gx, Gy, Gz);
                fidCellArray = [fidCellArray, fid];
            case 'acquirefov' % {'acquireFOV', res, FOV, acqTime, axis}
                FOV = curElement{3}; % mm
                res = curElement{2}; % Integer
                acqTime = curElement{4}; % ms
                dx = FOV/res;
                kMax = abs(1/dx); % 1/mm
                switch lower(curElement{5})
                    case 'x'
                        Gx = kMax/acqTime; % kHz/mm
                        Gy = 0; 
                        Gz = 0;
                    case 'y'
                        Gy = kMax/acqTime; % kHz/mm
                        Gz = 0; 
                        Gx = 0;
                    case 'z'
                        Gz = kMax/acqTime; % kHz/mm
                        Gx = 0; 
                        Gy = 0;
                    otherwise
                        error('Unrecognized axis %s in acquireFOV command', curElement{5});
                end
                [spins, fid] = AcquireGradRelax(spins, acqTime, res, Gx, Gy, Gz);
                fidCellArray = [fidCellArray, fid];
            otherwise
                fprintf('Unrecognized command in sequence: %s \n', cmd);
        end
    else % A pulse
        spins = ApplyPulseRelax(spins, curElement);
    end
end