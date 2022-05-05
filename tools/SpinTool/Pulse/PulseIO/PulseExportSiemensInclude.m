function pulseHeader = PulseExportSiemensInclude(pulse, filename, pulseBW, flipAngle, comment, minSlice, varName)
% SYNTAX: 
% 
%   pulseHeader = ExportPulseToIncludeFile(pulse, filename, ...
%                    pulseBW, flipAngle, comment, minSlice, varName)
% 
% Description: exports a given pulse structure (sans gradients) to
% a C++ include (.h) file, which can then be included in IDEA code
% and loaded into an arbitrary pulse structure. 
%
% Inputs:
%
% Variable Name   Units    Description
% pulse           -        Input RF pulse structure having the following
%                          fields:
%                          pulse.tp - duration of pulse, in ms
%                          pulse.RFamp - amplitude envelope (kHz)
%                          pulse.RFphase - phase of pulse as func. of time
%                          (radians)
%                          pulse.Gx - amplitude of x-gradient (in kHz/mm)
%                          pulse.Gy - amplitude of y-gradient (in kHz/mm)
%                          pulse.Gz - amplitude of z-gradient (in kHz/mm)
%                          Gx, Gy, Gz, RFamp, RFphase are all arrays with
%                          the same number of elements.
%                          If pulse is an array, such that pulse(k) is a
%                          pulse structure, this function will export a .h
%                          file with refGrad, etc ... being an array, and
%                          the amplitude & phase being 2D arrays. Furthermore,
%                          pulseBW, flipAngle both can be vectors now
%                          having the same number of elements as the pulse
%                          array (if left at a single number, it will be
%                          replicated for all pulses).
% filename        -        Output filename (including directory!) [1]
% pulseBW         kHz      Bandwidth of pulse, in kHz. Necessary for 
%                          calibrating reference gradient.
% flipAngle       deg.     Flip angle of pulse, in radians. Be sure to set
%                          the same flip angle in the Siemens sequence, using:
%                          myExternalPulse.setFlipAngle(flipAngle in degrees)
% [comment]       -        Comment to be added. Can be left blank (i.e.: ''),
%                          in which case 'Just a pulse' will be used. Same
%                          if omitted.
% [minSlice]      mm       Minimum excitable slice. If omitted, will be set
%                          to 1 mm [2]
% [varName]       -        String, used to define the variables in the
%                          include file. If omitted, the filename will be
%                          used (without the extension)[3] 
%
% [1] The filename will always have a .h extension (if omitted, or another
%     is provided, a .h will be placed automatically instead).
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
%
% [3] The .h file will have the pulse's attributes defined as static
%     variables. For example, suppose varName = 'myPulse'. Then
%     static double myPulseAmp = {0.0 0.1 0.4 .... };
%     static double myPulsePhase = {0.0 180.0 90.0 ...};  // in deg.
%     static double myPulseAmpInt = 42.1; 
%     static double myPulseAbsInt = ...
%     and so forth.

numPulses = numel(pulse);

if numel(pulseBW)==1
    pulseBW = ones(1, numPulses)*pulseBW;
else
    numel(pulseBW)
    if numel(pulseBW)~=numPulses
        error('Number of elements of pulseBW must either be 1 or equal to the number of pulses (%d)', numPulses);
    end
end

if numel(flipAngle)==1
    flipAngle = ones(1, numPulses)*flipAngle;
else
    if numel(flipAngle)~=numPulses
        error('Number of elements of flipAngle must either be 1 or equal to the number of pulses (%d)', numPulses);
    end
end

if numPulses>1
    for idxPulse=1:numPulses, numSteps(idxPulse) = numel(pulse(idxPulse).RFamp); end
    if any(diff(numSteps)), error('You have provided a pulse array with %d pulses. They all must have the same number of steps (they do not at the moment).', numPulses); end
else
    numSteps = numel(pulse.RFamp);
end

% Check filename extension. Remove if present. Append .h extension.
filenameNoExt = RemoveExtensionFromFilename(filename);
filename = [filenameNoExt,'.h'];

if nargin<6
    minSlice = 1;
end

if nargin<7
    varName = filenameNoExt;
end

% Check comment is not empty
if nargin<5
    comment = 'No user supplied comment.';
else
    if isempty(comment)
        comment = 'No user supplied comment.';
    end
end

% Open file
fid = fopen(filename,'w');

fprintf(fid,'// =------------------=\n');
fprintf(fid,'// Arbitrary pulse file\n');
fprintf(fid,'// =------------------=\n');
fprintf(fid,'//\n');
if numPulses>1
    fprintf(fid,'float %sRefGrad[%d];\n',  varName, numPulses);
    fprintf(fid,'float %sMinSlice[%d];\n', varName, numPulses);
    fprintf(fid,'float %sMaxSlice[%d];\n', varName, numPulses);
    fprintf(fid,'float %sAmpInt[%d];\n',   varName, numPulses);
    fprintf(fid,'float %sPowerInt[%d];\n', varName, numPulses);
    fprintf(fid,'float %sAbsInt[%d];\n',   varName, numPulses);
end

% Convert deg. --> rad.
flipAngle = flipAngle/180*pi;

for idxPulse=1:numPulses
    

    % Rescale input pulse. 
    % Magnitude is scaled to [0 1], Phase is scaled mod. 2 pi.
    maxAmp = max(abs(pulse(idxPulse).RFamp));  % in kHz
    if (maxAmp == 0)
        maxAmp = 1;
    end
    pulse(idxPulse).RFamp = pulse(idxPulse).RFamp./maxAmp;
    % If amplitude is negative for any reason, make it positive and add pi to
    % the phase
    pulse(idxPulse).RFphase = pulse(idxPulse).RFphase + ((1-sign(pulse(idxPulse).RFamp))/2)*pi;
    pulse(idxPulse).RFamp = abs(pulse(idxPulse).RFamp);
    % Make sure phaes is between [0, 2*pi], as demanded by Siemens
    pulse(idxPulse).RFphase = mod(pulse(idxPulse).RFphase, 2*pi);

    % Compute necessary integrals for Siemens .PTA file
    absIntegral = sum(pulse(idxPulse).RFamp);
    powerIntegral = sum(pulse(idxPulse).RFamp.^2);

    % Compute amplitude integral. See "External RF Pulses.doc" 
    % in the "IDEA - Educational\My Notes" folder for more details.
    dwellTime = pulse(idxPulse).tp/numSteps(idxPulse); % ms
    refAmp = 0.5;  % Amplitude of reference pulse is 0.5 kHz
    ampIntegral = (flipAngle(idxPulse)/pi) * (1/dwellTime) * (refAmp/maxAmp);

    % Compute the reference gradient.
    refGrad(idxPulse) = ComputeReferenceGradientSiemens(pulse(idxPulse).tp, pulseBW(idxPulse));

    % Create external file header
    fprintf(fid,'// =------------------=\n');
    fprintf(fid,'// Pulse %d / %d       \n', idxPulse, numPulses);
    fprintf(fid,'// =------------------=\n');
    fprintf(fid,'//\n');
    fprintf(fid,'// General statistics (as exported): \n');
    fprintf(fid,'//    Duration:     %.3f ms\n', pulse(idxPulse).tp);
    fprintf(fid,'//    Max B1:       %.3f kHz\n', max(pulse(idxPulse).RFamp));
    fprintf(fid,'//    Num. Steps:   %d\n', numSteps(idxPulse));
    fprintf(fid,'//    SAR*:         %.3f\n', CalcSAR(pulse(idxPulse), 'ref'));
    fprintf(fid,'//    Ref. Grad:    %.3f mT/m (based on %.3f kHz BW) \n', refGrad(idxPulse), pulseBW(idxPulse));
    fprintf(fid,'//    Ref. Grad:    %.3f kHz/mm \n', refGrad(idxPulse)*GetGyromagneticRatio('1h')/1000);
    fprintf(fid,'//    Flip Angle:   %.1f deg. \n', flipAngle(idxPulse)/pi*180);
    fprintf(fid,'// * - Relative to a 1 ms pi-pulse.\n');
    fprintf(fid,'//\n');
    fprintf(fid,'// User comment: %s\n', comment);
    fprintf(fid,'//\n');
    fprintf(fid,'// Hint: do not forget to include these global libraries and definitions: \n');
    fprintf(fid,'//    #include "MrServers\\MrMeasSrv\\SeqIF\\libRT\\libRT.h"\n');
    fprintf(fid,'//    #include "MrServers\\MrMeasSrv\\SeqFW\\libSSL\\SSL_local.h"\n');
    if numPulses>1
        fprintf(fid,'//    static sSample %sPulseArray[%d][%d];\n', varName, numPulses, numSteps(1));
    else
        fprintf(fid,'//    static sSample %sPulseArray[%d];\n', varName, numSteps);
    end
    fprintf(fid,'\n');
    if numPulses>1
        fprintf(fid,'%sRefGrad[%d]  = %0.4f;\n', varName, idxPulse-1, refGrad(idxPulse));
        fprintf(fid,'%sMinSlice[%d] = %.1f;\n', varName, idxPulse-1, minSlice);
        fprintf(fid,'%sMaxSlice[%d] = 200.0;\n', varName, idxPulse-1);
        fprintf(fid,'%sAmpInt[%d]   = %.4f;\n', varName,  idxPulse-1, ampIntegral);
        fprintf(fid,'%sPowerInt[%d] = %.4f;\n', varName, idxPulse-1, powerIntegral);
        fprintf(fid,'%sAbsInt[%d]   = %.4f;\n', varName, idxPulse-1, absIntegral);
    else
        fprintf(fid,'float %sRefGrad = %0.4f;\n', varName, refGrad(idxPulse));
        fprintf(fid,'float %sMinSlice = %.1f;\n', varName, minSlice);
        fprintf(fid,'float %sMaxSlice = 200.0;\n', varName);
        fprintf(fid,'float %sAmpInt = %.4f;\n', varName, ampIntegral);
        fprintf(fid,'float %sPowerInt = %.4f;\n', varName, powerIntegral);
        fprintf(fid,'float %sAbsInt = %.4f;\n', varName, absIntegral);
    end
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    for idx=1:numSteps
        if numPulses>1
            fprintf(fid,'%sPulseArray[%d][%d].flAbs = float(%.6f);    %sPulseArray[%d][%d].flPha = float(%.6f);\n', varName,  idxPulse-1, idx-1, pulse(idxPulse).RFamp(idx), varName, idxPulse-1, idx-1, pulse(idxPulse).RFphase(idx));
        else
            fprintf(fid,'%sPulseArray[%d].flAbs = float(%.6f);    %sPulseArray[%d].flPha = float(%.6f);\n', varName, idx-1, pulse(idxPulse).RFamp(idx), varName, idx-1, pulse(idxPulse).RFphase(idx));
        end
    end

    % Save data to output variable
    pulseHeader(idxPulse).comment = comment;
    pulseHeader(idxPulse).refGrad = refGrad(idxPulse);
    pulseHeader(idxPulse).refGradKHzmm = refGrad(idxPulse)*GetGyromagneticRatio('1h')/1000;
    pulseHeader(idxPulse).ampIntegral = ampIntegral;
    pulseHeader(idxPulse).powerIntegral = powerIntegral;
    pulseHeader(idxPulse).absIntegral = absIntegral;
    pulseHeader(idxPulse).maxGradSlewRate = CalcMaxGradSlewRate(pulse(idxPulse));
    pulseHeader(idxPulse).maxRFSlewRate = CalcMaxRFSlewRate(pulse(idxPulse));
    pulseHeader(idxPulse).maxRFSlewRate = CalcSAR(pulse(idxPulse));
end

% Close file.
fclose(fid);