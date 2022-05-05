function [pulse, pulseHeader] = PulseReadSiemensInclude(filename, pulseDuration, flipAngle, sliceThickness, gradientAxis, centerFreq, gradRiseTime)
% SYNTAX: 
%
%     [pulse, pulseHeader] = ReadSiemensIncludeFile(filename, [pulseDuration=1], 
%                                                   [flipAngle=90], [sliceThickness], [gradientAxis]) 
%
%
% Reads a pulse from a .h file.
%
% Inputs:
%
% Variable Name   Units    Description
% filename        -        Full filename, including directory. .PTA 
%                          extension not required.
% pulseDuration   ms       Duration of pulse. OPTIONAL (default: 1).
% flipAngle       deg.     Flip angle, in degrees. OPTIONAL (default: 90).
% variablePrefix  -        Optional. Prefix of variable names. If omitted,
%                          the filename will be used. For example, if 
%                          'variablePrefix' = 'myPulse', the function will
%                          look for myPulseRefGrad, myPulseMinSlice, etc.
%                          in the .h file
% sliceThickness  mm       Optional. Thickness of slice to be excited
% gradientAxis    -        Optional. 'x', 'y', or 'z' (case insensitive)
% centerFreq      kHz      Optional (default: 0). If provided, the center
%                          of the pulse will be shifted to the frequency
%                          in question.
% gradRiseTime    ms       Rise time of gradients. If omitted or set to
%                          0, no rise/fall times will be used and gradients
%                          will be rectangular.
%
% Outputs:
%
% Variable Name   Units    Description
% pulseHeader     -        Structure containing info from the first few
%                          lines of the .PTA file. The following fields
%                          are included:
%                          pulseName
%                          comment      
%                          refGrad       mT/m (for a 5.12 ms pulse over 10 mm slice)
%                          minSlice      mm
%                          maxSlice      mm
%                          ampInt        (a.u.)
%                          powerInt      (a.u.)
%                          absInt        (a.u.)
%                          duration      ms
%                          pulseBW       kHz
%                          For their meaning, consult the IDEA manual.
%                          Furthermore, two additional fields are provided:
%                          ampVec   - normalized vec. containing amplitudes
%                          phaseVec - vector containing phases
% pulse           -        Pulse structure. 


if nargin<8
    gradUnits = 'khzmm';
end

% If pulseDuration and flipAngle are not specified, set to default
if nargin<7
    gradRiseTime = 0;
end

if nargin<6
    centerFreq = 0;
end

if (nargin<5)
    gradientAxis = 'z';
end

if (nargin<4)
    sliceThickness = 10;
end

if (nargin<3) 
    flipAngle = 90;
end

if (nargin<2)
    pulseDuration = 1;
end

% Check for .h extension. Append if missing
filename = RemoveExtensionFromFilename(filename);
filename = [filename, '.h'];

if ~exist(filename, 'file'),
    error('Cannot find file %s', filename);
end

% -----------------------------------------------------------------------
% Load pulse data
% ----------------------------------------------------------------------

% Initialize header
pulseHeader.pulseName = 'MyPulse.Pulse';
pulseHeader.comment = 'No comment.';
pulseHeader.refGrad = 1;
pulseHeader.minSlice = 1;  % mm
pulseHeader.maxSlice = 100;  % mm
pulseHeader.ampInt = 50.000; 
pulseHeader.powerInt = 0.000;
pulseHeader.absInt = 0.000;

% Load data
fid=fopen(filename, 'r'); 
pulseFileStr = fread(fid, inf, 'uint8=>char')';
fclose(fid);

[pulseHeader.refGrad, ok] = ReadParameterFromASCII(pulseFileStr, 'RefGrad', '=', 1);
[pulseHeader.minSlice, ok] = ReadParameterFromASCII(pulseFileStr, 'MinSlice', '=', 1);
[pulseHeader.maxSlice, ok] = ReadParameterFromASCII(pulseFileStr, 'MaxSlice', '=', 1);
[pulseHeader.ampInt, ok] = ReadParameterFromASCII(pulseFileStr, 'AmpInt', '=', 1);
[pulseHeader.powerInt, ok] = ReadParameterFromASCII(pulseFileStr, 'PowerInt', '=', 1);
[pulseHeader.absInt, ok] = ReadParameterFromASCII(pulseFileStr, 'AbsInt', '=', 1);

if ~ok
    % This means that there's probably an array of pulses in the file.
    % We'll get the number of pulses by looking at the size of the C++
    % array RefGrad[XX]
    idxL = strfind(pulseFileStr, 'RefGrad[')+8;
    idxR = strfind(pulseFileStr, '];')-1;
    numPulses = str2double(pulseFileStr(idxL(1):idxR(1)));
else
    numPulses = 1;
end

if numel(flipAngle)==1
    flipAngle = ones(1, numPulses)*flipAngle;
end

if numel(flipAngle)~=numPulses
    error('Number of provided flip angles (%d) must equal 1 or match the number of pulses (%d)', numel(flipAngle), numPulses);
end

for idxPulse=1:numPulses

    if numPulses==1
        initStr = 'PulseArray[0].flAbs';
    else
        initStr = sprintf('PulseArray[%d][0].flAbs', idxPulse-1);
    end
    
    if numPulses>1
        [pulseHeader(idxPulse).refGrad, ok]  = ReadParameterFromASCII(pulseFileStr, sprintf('RefGrad[%d]',  idxPulse-1), '=', 1);
        [pulseHeader(idxPulse).minSlice, ok] = ReadParameterFromASCII(pulseFileStr, sprintf('MinSlice[%d]', idxPulse-1), '=', 1);
        [pulseHeader(idxPulse).maxSlice, ok] = ReadParameterFromASCII(pulseFileStr, sprintf('MaxSlice[%d]', idxPulse-1), '=', 1);
        [pulseHeader(idxPulse).ampInt, ok]   = ReadParameterFromASCII(pulseFileStr, sprintf('AmpInt[%d]',   idxPulse-1), '=', 1);
        [pulseHeader(idxPulse).powerInt, ok] = ReadParameterFromASCII(pulseFileStr, sprintf('PowerInt[%d]', idxPulse-1), '=', 1);
        [pulseHeader(idxPulse).absInt, ok]   = ReadParameterFromASCII(pulseFileStr, sprintf('AbsInt[%d]',   idxPulse-1), '=', 1);
    end
    
    % Some of the parameters may be surrounded by text. E.g., ampint is usually
    % written as ampint = float(32.41). We next strip the text.
    if ischar(pulseHeader(idxPulse).refGrad),  idxDigits = regexp(pulseHeader(idxPulse).refGrad, '[0-9]');  pulseHeader.refGrad  = str2num(pulseHeader(idxPulse).refGrad(idxDigits(1):idxDigits(end)));  end
    if ischar(pulseHeader(idxPulse).minSlice), idxDigits = regexp(pulseHeader(idxPulse).minSlice, '[0-9]'); pulseHeader.minSlice = str2num(pulseHeader(idxPulse).minSlice(idxDigits(1):idxDigits(end))); end
    if ischar(pulseHeader(idxPulse).maxSlice), idxDigits = regexp(pulseHeader(idxPulse).maxSlice, '[0-9]'); pulseHeader.maxSlice = str2num(pulseHeader(idxPulse).maxSlice(idxDigits(1):idxDigits(end))); end
    if ischar(pulseHeader(idxPulse).ampInt),   idxDigits = regexp(pulseHeader(idxPulse).ampInt, '[0-9]');   pulseHeader.ampInt   = str2num(pulseHeader(idxPulse).ampInt(idxDigits(1):idxDigits(end)));   end
    if ischar(pulseHeader(idxPulse).powerInt), idxDigits = regexp(pulseHeader(idxPulse).powerInt, '[0-9]'); pulseHeader.powerInt = str2num(pulseHeader(idxPulse).powerInt(idxDigits(1):idxDigits(end))); end
    if ischar(pulseHeader(idxPulse).absInt),   idxDigits = regexp(pulseHeader(idxPulse).absInt, '[0-9]');   pulseHeader.absInt   = str2num(pulseHeader(idxPulse).absInt(idxDigits(1):idxDigits(end)));   end

    if numPulses>1
        tempIndices1 = regexp(pulseFileStr, ['PulseArray\[', num2str(idxPulse-1), '\]\[[0-9]+\].flAbs'], 'start');
        matchStr = ['PulseArray\[', num2str(idxPulse-1), '\]\[[0-9]+\].flPha = float([0-9?.]+\)'];
        tempIndices2 = regexp(pulseFileStr, matchStr, 'end');
        firstIndex = tempIndices1(1);
        lastIndex = tempIndices2(end)+1;
    else
        firstIndex = strfind(pulseFileStr, initStr);
        lastIndex = numel(pulseFileStr);
    end
    pulseDataStr = pulseFileStr(firstIndex:lastIndex);
    indices = strfind(pulseDataStr, 'float(');
    indicesClose = strfind(pulseDataStr, ');');
    odds = indices(1:2:end); % The odd indices are for the "left column", i.e., amplitude
    evens = indices(2:2:end); % The even indices are for the "right column", i.e., phase
    numPoints = numel(odds);

    % Calculate how many digits each point has.
    % I am implicitly assuming all time points have the same number of digits,
    % so I'm just using the first one.
    numDigits = indicesClose(1)-indices(1)-1; 

    ampVec = zeros(1, numPoints);
    phaseVec = zeros(1, numPoints);

    for idx=1:numPoints
        ampVec(idx) = str2double(pulseDataStr(odds(idx)+6:odds(idx)+numDigits));
        phaseVec(idx) = str2double(pulseDataStr(evens(idx)+6:evens(idx)+numDigits));
    end

    pulseHeader(idxPulse).ampVec = ampVec./max(abs(ampVec));  % Normalize to [0,1]
    pulseHeader(idxPulse).phaseVec = phaseVec;

    % -----------------------------------------------------------------------
    % Create & Calibrate Pulse Structure
    % ----------------------------------------------------------------------

    pulse(idxPulse) = PulseCreateFromSiemensParameters(...
        pulseHeader(idxPulse).ampVec, ...
        pulseHeader(idxPulse).phaseVec, ...
        pulseHeader(idxPulse).ampInt, ...
        pulseHeader(idxPulse).refGrad, ...
        gradientAxis, ...
        pulseDuration, ...
        sliceThickness, ...
        flipAngle(idxPulse));

    % Add rise & fall times
    isSpatiallySelectivePulse = (sum(abs(pulse(idxPulse).Gx) + abs(pulse(idxPulse).Gy) + abs(pulse(idxPulse).Gz))>0);
    if gradRiseTime>0 && isSpatiallySelectivePulse
        pulse(idxPulse) = PulseAddGradRamp(pulse(idxPulse), gradientAxis, gradRiseTime);
    end


    pulseHeader(idxPulse).maxGradSlewRate = CalcMaxGradSlewRate(pulse(idxPulse));
    pulseHeader(idxPulse).maxRFSlewRate = CalcMaxRFSlewRate(pulse(idxPulse));

    refSliceThickness = 10; % mm
    pulseHeader(idxPulse).duration = pulse(idxPulse).tp; % ms
    refDuration = 5.12; % ms
    pulseHeader(idxPulse).pulseBW = pulseHeader(idxPulse).refGrad*GetGyromagneticRatio('1h')/1000*refSliceThickness*(refDuration/pulse(idxPulse).tp);

    % -----------------------------------------------------------------------
    % Shift pulse center in frequency space
    % ----------------------------------------------------------------------

    if centerFreq~=0
        pulse(idxPulse) = PulseShiftOffset(pulse(idxPulse), centerFreq); 
    end
end