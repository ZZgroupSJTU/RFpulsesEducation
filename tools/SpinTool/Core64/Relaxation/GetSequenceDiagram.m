function [tAxis, B1Amp, B1Phase, Gx, Gy, Gz, pulseCenters, pulseNames, GxSS, GySS, GzSS, GxPurge, GyPurge, GzPurge] = GetSequenceDiagram(seq, varargin)
% GetSequenceDiagram  Calculates the full sequence diagram
%   [tAxis, B1Amp, B1Phase, Gx, Gy, Gz] = GetSequenceDiagram(seq) Returns
%   B1(t), G(t) as a function of time (tAxis).
%
%   [tAxis, B1Amp, B1Phase, Gx, Gy, Gz] = GetSequenceDiagram
%   (seq, 'paramname', paramvalue) Allows the user to set several sequence
%   plotting options:
%
%     'HardPulseDuration'  Allows to override the default 0.1 us hard pulse
%                          duration for visualization purposes. This will 
%                          NOT alter the amplitude of the hard pulse, which 
%                          you can scale using HardPulseScaling. 
%                          Default is 0.1 us. Units are in ms!
%     'HardPulseScaling'   Scales hard pulse amplitudes by the given scaling
%                          factor. Default: 1 (no scaling). Alternatively, 
%                          you could use HardPulse90Amp.
%     'HardPulse90Amp'     An alternative to HardPulseScaling, this allows
%                          the user to set the amplitude of a hard 90-pulse 
%                          in kHz. Scales all other hard pulses accordingly. 
%                          Set to 0 just use whatever physical value is 
%                          calculated by the program. Default: 0.
%     'KillMxyDuration'    Default: 0. If set to 0, 'killmxy' events
%                          will not be drawn. Any other value signifies the 
%                          duration of the killmxy event, for visualization 
%                          purposes only! 
%     'AMPulsesZeroPhase'  Default: true. If set to true, the B1 amplitude
%                          of all amplitude modulated pulses (excluding
%                          hard pulses) will assume negative values, which
%                          represent 180-deg. phases. If set to false, 
%                          B1 amplitudes will be kept >=0 and negative
%                          amplitudes will be represented by a 180-deg.
%                          phase.
%     'PurgeRampTime'      Fixed rise time for gradient rise/fall times. 
%                          Set to 0 for no rise/fall times. Default: 0
%     'Accuracy'           Number of accuracy in digits when displaying
%                          pulse names, if degrees are used. Default: 0.

p = inputParser;


%
p.addParameter('HardPulseDuration', 0.0001, @(x) isa(x, 'double'));  
p.addParameter('HardPulseScaling', 1, @(x) isa(x, 'double'));
p.addParameter('HardPulse90Amp', 0, @(x) isa(x, 'double'));
p.addParameter('AMPulsesZeroPhase', true, @(x) isa(x, 'logical'));
p.addParameter('KillMxyDuration', 0, @(x) isa(x, 'double'));
p.addParameter('PurgeRampTime', 0, @(x) isa(x,'double'));
p.addParameter('Accuracy', 0, @(x) isnumeric(x) && (x>=0));
p.addParameter('ShowPulsePhase', true, @(x) islogical(x));
p.parse(varargin{:});

% 'HardPulse90Amp', 1, 'HardPulseDuration', 1, 'PurgeRampTime', 1)

freqADC = 0;

hardPulseScaling = p.Results.HardPulseScaling;
hardPulse90Amp = p.Results.HardPulse90Amp;
hardPulseDuration = p.Results.HardPulseDuration;
isAMPulsesZeroPhase = p.Results.AMPulsesZeroPhase;
purgeRampTime = p.Results.PurgeRampTime;
accuracy = p.Results.Accuracy;
isShowPulsePhase = p.Results.ShowPulsePhase;


pulseCenters = [];
pulseNames = {};

numElements = numel(seq);

% Initialize all time vectors
tAxis = 0;
B1Amp = 0;
B1Phase = 0;
Gx = 0;
Gy = 0;
Gz = 0;
GxSS = 0;
GySS = 0;
GzSS = 0;
GxPurge = 0;
GyPurge = 0;
GzPurge = 0;

% Calculate the maximal B1 amplitude (in kHz) during the sequence, including
% the hard pulses.
maxB1 = 0;
for idxEl=1:numElements
    curElement = seq{idxEl};
    if iscell(curElement)
        cmd = curElement{1};
        switch lower(cmd)
            case 'delay' % {'delay', d, [numSteps]}
            case 'hard' % {'hard', ang, ph} 
                if hardPulse90Amp~=0
                    maxB1 = max(maxB1, hardPulse90Amp);
                else
                    flipAngle = curElement{2};
                    pulseAmp90 = 1/(4*hardPulseDuration); % kHz
                    pulseAmp = flipAngle/90*pulseAmp90; 
                    maxB1 = max(maxB1, pulseAmp);
                end
            case 'rect' % {'rect', ang, ph, d, [numSteps]}   
                flipAngle = curElement{2};
                pulsePhase = curElement{3};
                duration = curElement{4};
                B1kHz = flipAngle/(360*duration);
                maxB1 = max(maxB1, B1kHz);
            case 'killmxy'  % {'killmxy'}
            case 'purge' % {'purge', Gx, Gy, Gz, d}
            case 'purgemoment' % {'purgemoment', kx, ky, kz, d}
            case 'acquire' % {'acquire', Nt, SW, Gx, Gy, Gz}
            case 'acquirefov' % {'acquireFOV', res, FOV, acqTime, axis}
            case 'pulse'
                curPulse = curElement{2};
                maxB1 = max(maxB1, max(curPulse.RFamp));
            otherwise
                fprintf('GetSequenceDiagram Warning: Unrecognized command in sequence: %s \n', cmd);
        end
    else % A pulse
        curPulse = curElement;
        maxB1 = max(maxB1, max(curPulse.RFamp));
    end
end


for idxEl=1:numElements
    curElement = seq{idxEl};
    if iscell(curElement)
        cmd = curElement{1};
        numVars = numel(curElement);
        switch lower(cmd)
            case 'delay' % {'delay', d, [numSteps]}
                if numel(curElement)<=2
                    NN = 1;
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
                flipAngle = curElement{2};
                if hardPulse90Amp>0
                    realHardPulse90Amp = 1/(4*hardPulseDuration); % kHz
                    flipAngleCalib = flipAngle*(hardPulse90Amp/realHardPulse90Amp);
                else
                    flipAngleCalib = flipAngle;
                end
                if numVars<3
                    pulsePhase = 0;
                else
                    pulsePhase = curElement{3};
                end
                numHardSteps = 20;
                curPulse = PulseCreateHard(hardPulseDuration*1000, flipAngleCalib*hardPulseScaling, pulsePhase, numHardSteps); % PulseCreateHard takes the duration in us!
            case 'rect' % {'rect', ang, ph, d, [numSteps]}   
                flipAngle = curElement{2};
                pulsePhase = curElement{3};
                if numVars<=4
                    curPulse = PulseCreateHard(curElement{4}*1000, flipAngle, pulsePhase);
                else
                    pulseDuration = curElement{4}; % ms
                    numSteps = curElement{5}; 
                    pulsePhase = curElement{3}/180*pi; % radians
                    peakB1 = flipAngle/360/pulseDuration; % kHz
                    curPulse = PulseCreateConst(pulseDuration, numSteps, peakB1, pulsePhase);
                end
            case 'killmxy'  % {'killmxy'}
                tAxis = [tAxis, tAxis(end)+1e-6];
                curPulse = [];
            case 'purge' % {'purge', Gx, Gy, Gz, d}
                GxAmp = curElement{2}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                GyAmp = curElement{3}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                GzAmp = curElement{4}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                flatTime = curElement{5};
                if purgeRampTime>0
                    NN = 50;
                    dtFlat = flatTime/NN;
                    numPurgeRampSteps = ceil(purgeRampTime/dtFlat);
                    newGradRampTime = numPurgeRampSteps*dtFlat;
                    dGx = GxAmp/numPurgeRampSteps;
                    dGy = GyAmp/numPurgeRampSteps;
                    dGz = GzAmp/numPurgeRampSteps;
                    curPulse.tp = flatTime + 2*newGradRampTime;
                    curPulse.RFamp   = zeros(1,NN + 2*numPurgeRampSteps); 
                    curPulse.RFphase = zeros(1,NN + 2*numPurgeRampSteps);
                    curPulse.Gx      = [0:dGx:GxAmp-dGx, GxAmp*ones(1,NN), GxAmp-dGx:-dGx:0];
                    curPulse.Gy      = [0:dGy:GyAmp-dGy, GyAmp*ones(1,NN), GyAmp-dGy:-dGy:0];
                    curPulse.Gz      = [0:dGz:GzAmp-dGz, GzAmp*ones(1,NN), GzAmp-dGz:-dGz:0];
                else
                    NN = 1;
                    curPulse.tp = flatTime;
                    curPulse.RFamp   = zeros(1,NN); 
                    curPulse.RFphase = zeros(1,NN);
                    curPulse.Gx      = GxAmp*ones(1,NN);
                    curPulse.Gy      = GyAmp*ones(1,NN);
                    curPulse.Gz      = GzAmp*ones(1,NN);
                end
            case 'purgemoment' % {'purgemoment', kx, ky, kz, d}
                kx = curElement{2}; % m^(-1)
                ky = curElement{3}; % m^(-1)
                kz = curElement{4}; % m^(-1)
                flatTime = curElement{5}; % ms
                GxAmp = 0.001*kx/flatTime; % kHz/mm
                GyAmp = 0.001*ky/flatTime; % kHz/mm
                GzAmp = 0.001*kz/flatTime; % kHz/mm
                if purgeRampTime>0
                    NN = 50;
                    dtFlat = flatTime/NN;
                    numPurgeRampSteps = ceil(purgeRampTime/dtFlat);
                    newGradRampTime = numPurgeRampSteps*dtFlat;
                    dGx = GxAmp/numPurgeRampSteps;
                    dGy = GyAmp/numPurgeRampSteps;
                    dGz = GzAmp/numPurgeRampSteps;
                    curPulse.tp = flatTime + 2*newGradRampTime;
                    curPulse.RFamp   = zeros(1,NN + 2*numPurgeRampSteps); 
                    curPulse.RFphase = zeros(1,NN + 2*numPurgeRampSteps);
                    curPulse.Gx      = [0:dGx:GxAmp-dGx, GxAmp*ones(1,NN), GxAmp-dGx:-dGx:0];
                    curPulse.Gy      = [0:dGy:GyAmp-dGy, GyAmp*ones(1,NN), GyAmp-dGy:-dGy:0];
                    curPulse.Gz      = [0:dGz:GzAmp-dGz, GzAmp*ones(1,NN), GzAmp-dGz:-dGz:0];
                else
                    NN = 1;
                    curPulse.tp = flatTime;
                    curPulse.RFamp   = zeros(1,NN); 
                    curPulse.RFphase = zeros(1,NN);
                    curPulse.Gx      = GxAmp*ones(1,NN);
                    curPulse.Gy      = GyAmp*ones(1,NN);
                    curPulse.Gz      = GzAmp*ones(1,NN);
                end
            case 'acquire' % {'acquire', Nt, SW, Gx, Gy, Gz}
                numAcqPoints = curElement{2};
                SW = curElement{3};
                dt = 1/SW;
                acqTime = dt*numAcqPoints;
                GxAmp = curElement{4}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                GyAmp = curElement{5}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                GzAmp = curElement{6}*GetGyromagneticRatio('1h')/1000; % Convert to kHz/mm
                tt = [0:dt:acqTime-dt];
                ttNorm = tt./acqTime;
                curPulse.tp = acqTime;
                curPulse.RFamp   = maxB1*cos(ttNorm*2*pi*freqADC).*exp(-ttNorm*5);
                curPulse.RFphase = zeros(1,numAcqPoints);
                curPulse.Gx      = GxAmp*ones(1,numAcqPoints);
                curPulse.Gy      = GyAmp*ones(1,numAcqPoints);
                curPulse.Gz      = GzAmp*ones(1,numAcqPoints);
            case 'acquirefov' % {'acquireFOV', res, FOV, acqTime, axis}
                FOV = curElement{3}; % mm
                res = curElement{2}; % Integer
                acqTime = curElement{4}; % ms
                dx = FOV/res;
                kMax = 1/dx; % 1/mm
                switch lower(curElement{5})
                    case 'x'
                        GxAmp = kMax/acqTime; % kHz/mm
                        GyAmp = 0; 
                        GzAmp = 0;
                    case 'y'
                        GyAmp = kMax/acqTime; % kHz/mm
                        GzAmp = 0; 
                        GxAmp = 0;
                    case 'z'
                        GzAmp = kMax/acqTime; % kHz/mm
                        GxAmp = 0; 
                        GyAmp = 0;
                    otherwise
                        error('Unrecognized axis %s in acquireFOV command', curElement{5});
                end
                dt = acqTime/res;
                tt = [0:dt:acqTime-dt];
                ttNorm = tt./acqTime;
                curPulse.tp = acqTime;
                curPulse.RFamp   = maxB1*cos(ttNorm*2*pi*freqADC).*exp(-ttNorm*5);
                curPulse.RFphase = zeros(1,res);
                curPulse.Gx      = GxAmp*ones(1,res);
                curPulse.Gy      = GyAmp*ones(1,res);
                curPulse.Gz      = GzAmp*ones(1,res);
            case 'pulse'
                curPulse = curElement{2};
                if numel(curElement)<3
                    curPhase = 0;
                else
                    curPhase = curElement{3};
                end
                if IsPulseAM(curPulse) && isAMPulsesZeroPhase
                    pulsePhase = mod(curPulse.RFphase, 2*pi);
                    errorThreshold = 1e-4;
                    curPulse.RFamp = -curPulse.RFamp.*((-1).^((pulsePhase)<errorThreshold));
                    curPulse.RFphase = curPulse.RFphase.*0 + curPhase;
                end
            otherwise
                fprintf('GetSequenceDiagram Warning: Unrecognized command in sequence: %s \n', cmd);
        end
    else % A pulse
        curPulse = curElement;
        if IsPulseAM(curPulse) && isAMPulsesZeroPhase
            pulsePhase = mod(curPulse.RFphase, 2*pi);
            errorThreshold = 1e-4;
            curPulse.RFamp = -curPulse.RFamp.*((-1).^((pulsePhase)<errorThreshold));
            curPulse.RFphase = curPulse.RFphase.*0;
        end
    end
    if ~isempty(curPulse)
        dt = curPulse.tp/numel(curPulse.RFamp);
        tAxisNew = tAxis(end) + [dt:dt:curPulse.tp];
        if isempty(tAxis)
            tAxis = [0, tAxisNew, tAxisNew(end), tAxisNew(end)]; 
        else
            tAxis = [tAxis, tAxis(end), tAxisNew, tAxisNew(end), tAxisNew(end)]; 
        end
        B1Amp = [B1Amp, 0, curPulse.RFamp, curPulse.RFamp(end), 0];
        B1Phase = [B1Phase, 0, curPulse.RFphase, curPulse.RFphase(end), 0];
        Gx = [Gx, 0, curPulse.Gx, curPulse.Gx(end), 0];
        Gy = [Gy, 0, curPulse.Gy, curPulse.Gy(end), 0];
        Gz = [Gz, 0, curPulse.Gz, curPulse.Gz(end), 0];
        if iscell(curElement)
            switch lower(cmd)
                case {'purge', 'purgemoment'}
                    GxSS = [GxSS, zeros(1,numel(curPulse.Gx)+3)];
                    GySS = [GySS, zeros(1,numel(curPulse.Gy)+3)];
                    GzSS = [GzSS, zeros(1,numel(curPulse.Gz)+3)];
                    GxPurge = [GxPurge, 0, curPulse.Gx, curPulse.Gx(end), 0];
                    GyPurge = [GyPurge, 0, curPulse.Gy, curPulse.Gy(end), 0];
                    GzPurge = [GzPurge, 0, curPulse.Gz, curPulse.Gz(end), 0];
                case {'pulse'}
                    GxSS = [GxSS, 0, curPulse.Gx, curPulse.Gx(end), 0];
                    GySS = [GySS, 0, curPulse.Gy, curPulse.Gy(end), 0];
                    GzSS = [GzSS, 0, curPulse.Gz, curPulse.Gz(end), 0];
                    GxPurge = [GxPurge, zeros(1,numel(curPulse.Gx)+3)];
                    GyPurge = [GyPurge, zeros(1,numel(curPulse.Gy)+3)];
                    GzPurge = [GzPurge, zeros(1,numel(curPulse.Gz)+3)];
                otherwise
                    GxSS = [GxSS, zeros(1,numel(curPulse.Gx)+3)];
                    GySS = [GySS, zeros(1,numel(curPulse.Gy)+3)];
                    GzSS = [GzSS, zeros(1,numel(curPulse.Gz)+3)];
                    GxPurge = [GxPurge, zeros(1,numel(curPulse.Gx)+3)];
                    GyPurge = [GyPurge, zeros(1,numel(curPulse.Gy)+3)];
                    GzPurge = [GzPurge, zeros(1,numel(curPulse.Gz)+3)];
            end
        else
            GxSS = [GxSS, 0,curPulse.Gx, curPulse.Gx(end), 0];
            GySS = [GySS, 0,curPulse.Gy, curPulse.Gy(end), 0];
            GzSS = [GzSS, 0,curPulse.Gz, curPulse.Gz(end), 0];
            GxPurge = [GxPurge, zeros(1,numel(curPulse.Gx)+3)];
            GyPurge = [GyPurge, zeros(1,numel(curPulse.Gy)+3)];
            GzPurge = [GzPurge, zeros(1,numel(curPulse.Gz)+3)];
        end
        
        if iscell(curElement)
            switch lower(cmd)
                case {'hard', 'rect'}
                    pulseCenters = [pulseCenters, (tAxisNew(end) + tAxisNew(1))/2];
                    switch pulsePhase
                        case 0
                            phaseStr = 'x';
                        case 90
                            phaseStr = 'y';
                        case 180
                            phaseStr = '-x';
                        case 270
                            phaseStr = '-y';
                        otherwise
                            phaseStr = sprintf(['%.', num2str(accuracy),'f°'], pulsePhase);
                    end
                    if isShowPulsePhase
                        pulseNames{end+1} = sprintf(['%.', num2str(accuracy),'f°_{%s}'], flipAngle, phaseStr);
                    else
                        pulseNames{end+1} = sprintf(['%.', num2str(accuracy),'f°'], flipAngle);
                    end
            end
        else
            pulseCenters = [pulseCenters, (tAxisNew(end) + tAxisNew(1))/2];
            if IsPulseAM(curElement)
                flipAngle = curElement.tp/numel(curElement.RFamp)*abs(sum(curElement.RFamp.*exp(1i*curElement.RFphase)))*360; % In deg.
                flipAngle = round(flipAngle);
                pulsePhase = curElement.RFphase(1)/pi*180; % in Deg.
                pulsePhase = round(pulsePhase*10)/10; % First decimal point
                pulsePhase = mod(pulsePhase, 360);
                switch pulsePhase
                    case 0
                        phaseStr = 'x';
                    case 90
                        phaseStr = 'y';
                    case 180
                        phaseStr = '-x';
                    case 270
                        phaseStr = '-y';
                    otherwise
                        phaseStr = sprintf(['%.', num2str(accuracy),'f°'], pulsePhase);
                end
                if isShowPulsePhase
                    pulseNames{end+1} = sprintf(['%.', num2str(accuracy),'f°_{%s}'], flipAngle, phaseStr);
                else
                    pulseNames{end+1} = sprintf(['%.', num2str(accuracy),'f°'], flipAngle);
                end
            else
                pulseNames{end+1} = 'Pulse';
            end
        end
                
    end
end

B1Amp = [B1Amp, 0];
B1Phase = [B1Phase, 0];
Gx = [Gx, 0];
Gy = [Gy, 0];
Gz = [Gz, 0];
GxSS = [GxSS, 0];
GySS = [GySS, 0];
GzSS = [GzSS, 0];
GxPurge = [GxPurge, 0];
GyPurge = [GyPurge, 0];
GzPurge = [GzPurge, 0];
tAxis = [tAxis, 0];
