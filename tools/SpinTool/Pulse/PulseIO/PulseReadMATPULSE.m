function [pulse, pulseHeader] = PulseReadMATPULSE(filename, pulseDuration, maxB1)
% Reads a MATPULSE output text file
%   [pulse, pulseHeader] = PulseReadMATPULSE(filename, pulseDuration, maxB1)
%   Returns a pulse structure and some header information from the given
%   file (which is saved from the MATPULSE software). pulseDuration is in
%   ms and maxB1 is in kHz.

pulse = [];
pulseHeader = [];

if numel(maxB1)>1, error('maxB1 has more than one element.'); end
if numel(pulseDuration)>1, error('pulseDuration has more than one element.'); end

if ~exist(filename, 'file')
    error('Cannot find file %s', filename);
end


fid = fopen(filename);
counter = 0;
pulseAmp = [];
pulsePhase = [];
while true
    curLine = fgetl(fid);
    if ~ischar(curLine), break; end
    if strcmpi(curLine(1:2), '##')
        curLine = curLine(3:end);
        curData = regexp(curLine, '=', 'split');
        if numel(curData)==2
            curData = strtrim(curData);
            paramName = curData{1};
            paramName = regexprep(paramName, '[^\w]', ''); % Remove non-letter characters
            paramValue = curData{2};
            paramValueNum = str2double(paramValue);
            if ~isempty(paramValueNum) && ~isnan(paramValueNum)
                pulseHeader.(paramName) = paramValueNum;
            else
                pulseHeader.(paramName) = paramValue;
            end
        end
    else
        curData = str2num(curLine);
        if ~isempty(curData) && numel(curData)==2
            counter = counter + 1;
            pulseAmp(counter) = curData(1)/100*maxB1;
            pulsePhase(counter) = curData(2)/180*pi;
        end
    end
end
fclose(fid);
if counter==0, return; end

pulse.tp = pulseDuration;
pulse.RFamp = pulseAmp;
pulse.RFphase = pulsePhase;
pulse.Gx = zeros(1, counter);
pulse.Gy = zeros(1, counter);
pulse.Gz = zeros(1, counter);

