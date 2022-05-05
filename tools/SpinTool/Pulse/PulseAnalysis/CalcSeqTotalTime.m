function T = CalcSeqTotalTime(seq)
% T = CalcSeqTotalTime(seq)
%   Returns the total time of the given sequence.
%   For more information on sequences, consult ApplySequence.m
%
% Input Variables
% Variable Name   Size   Units   Default  Description
% seq             N/A    N/A     N/A      Input sequence. See ApplySequence
%                                         for more information.
% Output Variables
% Variable Name   Size   Units   Default  Description
% T               1x1    ms      N/A      Total sequence time

T = 0;
numElements = numel(seq);
for idxEl=1:numElements
    curElement = seq{idxEl};
    if iscell(curElement)
        cmd = curElement{1};
        switch lower(cmd)
            case 'pulse'
                pulse = curElement{2};
                T = T + pulse.tp;
            case 'delay' % {'delay', d}
                T = T + curElement{2};
            case 'hard' % {'hard', ang, ph}      
                T = T + 0.0001;
                % Do nothing
            case 'rect' % {'rect', ang, ph, d}   
                T = T + curElement{4};
            case 'killmxy'  % {'killmxy'}
                % Do nothing
            case 'purge' % {'purge', Gx, Gy, Gz, d}
                T = T + curElement{5};
            case 'purgemoment' % {'purgemoment', kx, ky, kz, d}
                T = T + curElement{5};
            case 'acquire' % {'acquire', Nt, SW, Gx, Gy, Gz}
                if ispulse(curElement{2})
                    T = T + curElement{2}.tp;
                else
                    numAcqPoints = curElement{2};
                    SW = curElement{3};
                    dt = 1/SW;
                    acqTime = dt*numAcqPoints;
                    T = T + acqTime;
                end
            case 'acquirefov' % {'acquireFOV', res, FOV, acqTime, axis}
                acqTime = curElement{4}; % ms
                T = T + acqTime; 
            case 'acquirebw' % {'acquireBW', res, FOV, BWPerPixel, axis}
                FOV = curElement{3}; % mm
                res = curElement{2}; % Integer
                BWPerPixel = curElement{4}*0.001; % Hz*0.001=kHz
                dx = FOV/res; % mm
                kMax = 1/dx; % 1/mm
                G = abs(BWPerPixel/dx); % kHz/mm 
                acqTime = abs(kMax/G); % ms
                T = T + acqTime; 
            otherwise
                fprintf('Unrecognized command in sequence: %s \n', cmd);
        end
    else % A pulse
        T = T + curElement.tp; 
    end
end
