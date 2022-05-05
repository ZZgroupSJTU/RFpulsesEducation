function p = PulseReadBinary(filename, duration, peakB1)
% SYNTAX:
%
%     p = ReadPulseBinary(filename, duration, peakB1)
%
% Reads a binary file containing a complex vector of numbers describing
% a pulse's x and y components, and returns a pulse structure
%
% Input Variables
% Variable Name    Units    Description
% filename         -        Full filename of binary file.
% [duration]       ms       Duration of pulse. If omitted, set to 1 ms.
% [peakB1]         kHz      Peak power of pulse. If omitted, the peak
%                           value as read for the file will be used instead.
%
% Output Variables
% Variable Name    Units    Description
% p                -        The pulse structure.
%
% Example: 
%
% % Reads the pulse from MyFile.bin and assigns it a duration of 10 ms
% % and a peak power of 2.5 kHz
% p = ReadPulseBinary('MyFile.bin', 10, 2.5);
%
%

if nargin<2, duration = 1; end

fid = fopen(filename, 'r');
pulseArray = fread(fid, inf, 'float32');
fclose(fid);

pulseX = pulseArray(1:2:end)';
pulseY = pulseArray(2:2:end)';
pulseComplex = pulseX + 1i*pulseY;
N = numel(pulseComplex);
pulseAmp = abs(pulseComplex);
pulsePhase = phase(pulseComplex);

if nargin>=3, 
    pulseAmp = pulseAmp./max(abs(pulseAmp)).*peakB1;
end

p.tp = duration;
p.RFamp = pulseAmp;
p.RFphase = pulsePhase;
p.Gx = zeros(1,N);
p.Gy = zeros(1,N);
p.Gz = zeros(1,N);


