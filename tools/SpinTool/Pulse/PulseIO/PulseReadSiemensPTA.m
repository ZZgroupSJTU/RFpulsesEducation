function [pulseHeader, pulse] = PulseReadSiemensPTA(filename, pulseDuration, flipAngle, sliceThickness, gradientAxis)
% Reads a pulse from a .PTA file.
%
% Inputs:
%
% Variable Name   Units    Description
% filename        -        Full filename, including directory. .PTA 
%                          extension not required.
% pulseDuration   ms       Duration of pulse. OPTIONAL (default: 1).
% flipAngle       deg.     Flip angle, in degrees. OPTIONAL (default: 90).
% sliceThickness  mm       Thickness of slice to be excited
% gradientAxis    -        'x', 'y', or 'z' (case insensitive)
%
%
% Outputs:
%
% Variable Name   Units    Description
% pulseHeader     -        Structure containing info from the first few
%                          lines of the .PTA file. The following fields
%                          are included:
%                          pulseName
%                          comment
%                          refGrad
%                          minSlice
%                          maxSlice
%                          ampInt
%                          powerInt
%                          absInt
%                          For their meaning, consult the IDEA manual.
%                          Furthermore, two additional fields are provided:
%                          ampVec   - normalized vec. containing amplitudes
%                          phaseVec - vector containing phases
% pulse           -        Pulse structure. 


% If pulseDuration and flipAngle are not specified, set to default
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

% Check for .PTA extension. Append if missing
filename = RemoveExtensionFromFilename(filename);
filename = [filename, '.pta'];


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
fid = fopen(filename, 'r'); 
phaseVec = [];
ampVec = [];
while (~feof(fid))
    curLine = fgetl(fid);
    % Strip comments
    commentBeginIndex = regexp(curLine, ';');
    if (~isempty(commentBeginIndex))
        curLine = curLine(1:commentBeginIndex-1);
    end
    % Check string is not empty
    if (~isempty(curLine))
        if (~isempty(str2num(curLine)))  % Try to convert current line to #
            numVec = str2num(curLine);
            ampVec = [ampVec, numVec(1)];
            phaseVec = [phaseVec, numVec(2)];
        else  % If current line doesn't contain only numbers, extract header info
            % Pulse name
            nameIdx = regexpi(curLine,'pulsename:');
            if (~isempty(nameIdx))
                curLine = strtrim(curLine(11:end));
                pulseHeader.pulseName = curLine;
            end
            % Comment
            commentIdx = regexpi(curLine,'comment:');
            if (~isempty(commentIdx))
                curLine = strtrim(curLine(9:end));
                pulseHeader.comment = curLine;
            end
            % Reference Gradient
            refGradIdx = regexpi(curLine,'refgrad:');
            if (~isempty(refGradIdx))
                curLine = strtrim(curLine(9:end));
                pulseHeader.refGrad = str2double(curLine);
            end
            % Minimum slice
            minSliceIdx = regexpi(curLine,'minslice:');
            if (~isempty(minSliceIdx))
                curLine = strtrim(curLine(10:end));
                pulseHeader.minSlice = str2double(curLine);
            end
            % Maximum slice
            maxSliceIdx = regexpi(curLine,'maxslice:');
            if (~isempty(maxSliceIdx))
                curLine = strtrim(curLine(10:end));
                pulseHeader.maxSlice = str2double(curLine);
            end
            % Amplitude Integral
            ampIntIdx = regexpi(curLine,'ampint:');
            if (~isempty(ampIntIdx))
                curLine = strtrim(curLine(8:end));
                pulseHeader.ampInt = str2double(curLine);
            end
            % Power Integral
            powerIntIdx = regexpi(curLine,'powerint:');
            if (~isempty(powerIntIdx))
                curLine = strtrim(curLine(10:end));
                pulseHeader.powerInt = str2double(curLine);
            end
            % Absolute Integral
            absIntIdx = regexpi(curLine,'absint:');
            if (~isempty(absIntIdx))
                curLine = strtrim(curLine(8:end));
                pulseHeader.absInt = str2double(curLine);
            end
            % Calculate (infer) bandwidth, in kHz
            pulseHeader.BW = pulseHeader.refGrad * 0.4257 * 5.12/pulseDuration; 
        end
    end
end


pulseHeader.ampVec = ampVec./max(abs(ampVec));  % Normalize to [0,1]
pulseHeader.phaseVec = phaseVec;

% -----------------------------------------------------------------------
% Create & Calibrate Pulse Structure
% ----------------------------------------------------------------------

pulse = PulseCreateFromSiemensParameters(...
    pulseHeader.ampVec, ...
    pulseHeader.phaseVec, ...
    pulseHeader.ampInt, ...
    pulseHeader.refGrad, ...
    gradientAxis, ...
    pulseDuration, ...
    sliceThickness, ...
    flipAngle);

pulseHeader.maxGradSlewRate = CalcMaxGradSlewRate(pulse);
pulseHeader.maxRFSlewRate = CalcMaxRFSlewRate(pulse);

% -----------------------------------------------------------------------
% Cleanup
% ----------------------------------------------------------------------

fclose(fid);