function [Mx, My, Mz, tAxis, B1Amp, B1Phase] = ApplySequenceDiagnostics(seq, numRep, cs, r, Min, T1, T2, M0)
% SYNTAX:
%
%    [Mx, My, Mz, tAxis, B1Amp, B1Phase] 
%           = ApplySequenceDiagnostics(seq, numRep, cs, r, Min, T1, T2, M0)
%
% Calculates Mx(t), My(t), Mz(t) as a function of time for a given sequence
% for a spin with well defined properties.
%
% Input Variables
% Variable Name   Size   Units   Default  Description
% seq             N/A    N/A     N/A      Input sequence. See ApplySequence
%                                         for more information.
% [numRep]        1x1    #       1        Number of repetitions of sequence
%                                         (e.g. to drive system into 
%                                         dynamic equilibrium) before 
%                                         starting to "acquire". If
%                                         positive, the repetitions will be
%                                         saved and shown. If negative, the
%                                         repetitions will be discarded.
% [cs]            1x1    kHz     0        Chemical shift
% [r]             3x1    mm      [0;0;0]  Position (x,y,z)
% [Min]           3x1    a.u.    [0;0;1]  Input magnetization
% [T1]            1x1    ms      1000     Longitudinal relaxation
% [T2]            1x1    ms      100      Transverse relaxation
% [M0]            1x1    a.u.    1        Equilibrium magnetization
%
% Output Variables
% Variable Name   Size   Units   Default  Description
% Mx, Mx, Mz      1xN    a.u.    N/A      Mx, My and Mz as a function of
%                                         time (i.e. as tAxis)
% tAxis           1xN    ms      N/A      Time axis
% B1Amp           1xN    kHz     N/A      Amplitude of B1(t) as a function
%                                         of tAxis.
% B1Phase         1xN    rad     N/A      Phase of B1(t) as a function of 
%                                         tAxis.

TR = CalcSeqTotalTime(seq);

if numRep==0 
    fprintf('Error in ApplySequenceDiagnostics: numRep = 0.\n');
    return
end

isDummyScans = (numRep<0);
numRep = abs(numRep);

% Concatenate pre-scans
if numRep>1
    seq = repmat(seq, 1, numRep);
end

numElements = numel(seq);

% Initialize time axis
tAxis = [0];

% Initialize magnetization
Mx = [Min(1)];
My = [Min(2)];
Mz = [Min(3)];

B1Amp = [0];
B1Phase = [0];

for idxEl=1:numElements
    curElement = seq{idxEl};
    if iscell(curElement)
        cmd = curElement{1};
        switch lower(cmd)
            case 'delay' % {'delay', d, [numSteps]}
                if numel(curElement)<=2
                    dt = min(T2,T1)/100;
                    NN = ceil(curElement{2}/dt);
                else
                    NN = curElement{3};
                end
                curPulse.tp = curElement{2};
                curPulse.RFamp   = zeros(1,NN); 
                curPulse.RFphase = zeros(1,NN);
                curPulse.Gx      = zeros(1,NN);
                curPulse.Gy      = zeros(1,NN);
                curPulse.Gz      = zeros(1,NN);
            case 'hard' % {'hard', ang, ph}
                if numel(curElement)<3
                    pulsePhase = 0;
                else
                    pulsePhase = curElement{3};
                end
                curPulse = PulseCreateHard(0.1, curElement{2}, pulsePhase);
            case 'rect' % {'rect', ang, ph, d, [numSteps]}   
                if numel(curElement)<=4
                    curPulse = PulseCreateHard(curElement{4}*1000, curElement{2}, curElement{3});
                else
                    pulseDuration = curElement{4}; % ms
                    numSteps = curElement{5}; 
                    pulsePhase = curElement{3}/180*pi; % radians
                    flipAngle = curElement{2}/360; % fraction of 360
                    peakB1 = flipAngle/pulseDuration; % kHz
                    curPulse = PulseCreateConst(pulseDuration, numSteps, peakB1, pulsePhase);
                end
            case 'killmxy'  % {'killmxy'}
                Mx = [Mx, 0];
                My = [My, 0];
                Mz = [Mz, Mz(end)];
                tAxis = [tAxis, tAxis(end)+1e-6];
                curPulse = [];
            case 'purge' % {'purge', Gx, Gy, Gz, d}
                Gx = curElement{2}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                Gy = curElement{3}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                Gz = curElement{4}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                tp = curElement{5};
                curPulse.tp = tp;
                dt = min(T2,T1)/100;
                NN = ceil(tp/dt);
                curPulse.RFamp   = zeros(1,NN); 
                curPulse.RFphase = zeros(1,NN);
                curPulse.Gx      = Gx*ones(1,NN);
                curPulse.Gy      = Gy*ones(1,NN);
                curPulse.Gz      = Gz*ones(1,NN);
            case 'purgemoment' % {'purgemoment', kx, ky, kz, d}
                kx = curElement{2}; % m^(-1)
                ky = curElement{3}; % m^(-1)
                kz = curElement{4}; % m^(-1)
                tp = curElement{5}; % ms
                Gx = 0.001*kx/tp; % kHz/mm
                Gy = 0.001*ky/tp; % kHz/mm
                Gz = 0.001*kz/tp; % kHz/mm
                curPulse.tp = tp;
                dt = min(T2,T1)/100;
                NN = ceil(tp/dt);
                curPulse.RFamp   = zeros(1,NN); 
                curPulse.RFphase = zeros(1,NN);
                curPulse.Gx      = Gx*ones(1,NN);
                curPulse.Gy      = Gy*ones(1,NN);
                curPulse.Gz      = Gz*ones(1,NN);
            case 'acquire' % {'acquire', Nt, SW, Gx, Gy, Gz}
                numAcqPoints = curElement{2};
                SW = curElement{3};
                dt = 1/SW;
                acqTime = dt*numAcqPoints;
                Gx = curElement{4}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                Gy = curElement{5}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                Gz = curElement{6}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                curPulse.tp = acqTime;
                curPulse.RFamp   = zeros(1,numAcqPoints); 
                curPulse.RFphase = zeros(1,numAcqPoints);
                curPulse.Gx      = Gx*ones(1,numAcqPoints);
                curPulse.Gy      = Gy*ones(1,numAcqPoints);
                curPulse.Gz      = Gz*ones(1,numAcqPoints);
            case 'pulse' 
                numVars = numel(curElement);
                curPulse = curElement{2};
                if numVars>2
                    curPulse.RFphase = curPulse.RFphase + curElement{3};
                end
            otherwise
                fprintf('Unrecognized command in sequence: %s \n', cmd);
        end
    else % A pulse
        curPulse = curElement;
    end
    if ~isempty(curPulse)
        Min = [Mx(end); My(end); Mz(end)];
        [MxTemp, MyTemp, MzTemp] = ApplyPulseRelaxDiagnostics(cs, r, Min, curPulse, T1, T2, M0);
        dt = curPulse.tp/numel(curPulse.RFamp);
        tAxis = [tAxis, tAxis(end)+[dt:dt:curPulse.tp]]; % We already know the magnetization at time t=0
        Mx = [Mx, MxTemp(2:end)];
        My = [My, MyTemp(2:end)];
        Mz = [Mz, MzTemp(2:end)];
        B1Amp = [B1Amp, curPulse.RFamp];
        B1Phase = [B1Phase, curPulse.RFphase];
    end
end

if isDummyScans
    % Discard "pre-scans"
    idx = find(tAxis>=TR*(numRep-1), 1, 'first');
    tMin = tAxis(idx);
    tAxis = tAxis(idx:end);
    tAxis = tAxis - tMin;
    Mx = Mx(idx:end);
    My = My(idx:end);
    Mz = Mz(idx:end);
end