function pulseHeader = PulseExportSiemensPTA(pulse, pulsename, filename, pulseBW, flipAngle, comment, minSlice)
% Description: exports a given pulse structure (sans gradients) to
% a siemens ASCII external file format. The Siemens format is detailed
% in section C.5-25 in the VB13 IDEA programming manual.
%
% Inputs:
%
% Variable Name   Units    Description
% pulse           -        Input RF pulse structure
% pulsename       -        Name of RF pulse
% filename        -        Output filename (including directory!) [1]
% pulseBW         kHz      Bandwidth of pulse, in kHz. Necessary for 
%                          calibrating reference gradient.
% flipAngle       rad.     Flip angle of pulse, in radians. Be sure to set
%                          the same flip angle in the Siemens sequence, using:
%                          myExternalPulse.setFlipAngle(flipAngle in degrees)
% comment         -        Comment to be added. Can be left blank (i.e.: ''),
%                          in which case 'Just a pulse' will be used.
% minSlice        mm       Minimum excitable slice. [2]
%
% [1] The filename must not include an extension. Any extension will be
% removed. A .PTA extension will be added (as demanded by Siemens).
%
% [2] A seemingly "harmless" quantity, but one that should be set by the
%     designer for the following reason: the pulse has a certain stepsize
%     (dwell time), which determines the overall affected bandwidth:
%     FOVBW = 1/(dwell time)
%     Beyond that, the spectral pattern of the pulse will replicate
%     itself. Thus, a pulse designed to excite a 1 cm slice with a 10 cm
%     FOV will excite a 1 mm slice with a 1 cm FOV if the gradient is 
%     multiplied by 10. As a result, if the object is 10 cm long, 
%     artifacts will arise by going below a 1 cm slice.

% Check filename extension. Remove if present. Append .PTA extension.
n = findstr(filename,'.');
if ~isempty(n)
    filename=filename(1:n-1);
end
filename = [filename,'.PTA'];

% Check comment is not empty
if isempty(comment)
    comment = 'Just a pulse.';
end

% Check angle is given in radians, not degrees!
if (flipAngle>2*pi)
    fprintf('POSSIBLE ERROR IN ExportPulseToSiemens\n');
    fprintf('You have entered flipAngle = %.2f in ExportPulseToSiemens. This routine expects an angle in radians.\n', flipAngle);
    flipAngle = flipAngle/360*2*pi;
    fprintf('Did you use degrees? I have changed your input to %.2f radians for you.\n', flipAngle);
end

% Number of time steps of RF pulse
numSteps = length(pulse.RFamp);

% Rescale input pulse. 
% Magnitude is scaled to [0 1], Phase is scaled mod. 2 pi.
maxAmp = max(abs(pulse.RFamp));  % in kHz
if (maxAmp == 0)
    maxAmp = 1;
end
pulse.RFamp = pulse.RFamp./maxAmp;
% If amplitude is negative for any reason, make it positive and add pi to
% the phase
pulse.RFphase = pulse.RFphase + ((1-sign(pulse.RFamp))/2)*pi;
pulse.RFamp = abs(pulse.RFamp);
% Make sure phaes is between [0, 2*pi], as demanded by Siemens
pulse.RFphase = mod(pulse.RFphase, 2*pi);

% Compute necessary integrals for Siemens .PTA file
absIntegral = sum(pulse.RFamp);
powerIntegral = sum(pulse.RFamp.^2);

% Compute amplitude integral. See "External RF Pulses.doc" 
% in the "IDEA - Educational\My Notes" folder for more details.
dwellTime = pulse.tp/numSteps;
refAmp = 0.5;  % Amplitude of reference pulse is 0.5 kHz
% ampIntegral = flipAngle / (2*pi * dwellTime * maxAmp);
ampIntegral = (flipAngle/pi) * (1/dwellTime) * (refAmp/maxAmp);

% Compute the reference gradient.
refGrad = ComputeReferenceGradientSiemens(pulse.tp, pulseBW);

% Open file
fid = fopen(filename,'w');

% Create external file header
fprintf(fid,'PULSENAME: %s \n',pulsename);
fprintf(fid,['COMMENT: ',comment,'\n']);
fprintf(fid,'REFGRAD: %0.8f\n',refGrad);
fprintf(fid,'MINSLICE: %d\n', minSlice); 
fprintf(fid,'MAXSLICE: 200.0\n');
fprintf(fid,'AMPINT: %0.8f\n', ampIntegral);
fprintf(fid,'POWERINT: %0.8f\n', powerIntegral);
fprintf(fid,'ABSINT: %0.8f\n', absIntegral);

% Save data to output variable
pulseHeader.pulsename = pulsename;
pulseHeader.comment = comment;
pulseHeader.refGrad = refGrad;
pulseHeader.refGradKHzmm = refGrad*GetGyromagneticRatio('1h')/1000;
pulseHeader.ampIntegral = ampIntegral;
pulseHeader.powerIntegral = powerIntegral;
pulseHeader.absIntegral = absIntegral;
pulseHeader.maxGradSlewRate = CalcMaxGradSlewRate(pulse);
pulseHeader.maxRFSlewRate = CalcMaxRFSlewRate(pulse);
pulseHeader.maxRFSlewRate = CalcSAR(pulse);

% Output RF data.
for idx=1:numSteps
    fprintf(fid,'%0.6f  %0.6f ; (%d)\n',pulse.RFamp(idx), pulse.RFphase(idx), idx);
end

% Close file.
fclose(fid);