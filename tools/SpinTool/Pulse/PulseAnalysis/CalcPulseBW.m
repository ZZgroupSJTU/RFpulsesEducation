function [BWBase, BWFWHM] = CalcPulseBW(pulse, threshold, pulseType, profileType)
% SYNTAX:
% 
%   [BWBase, BWFWHM] = CalcPulseBW(pulse, threshold, pulseType, profileType)
%
% Calculates the bandwidth of the pulse automatically, both FWHM and at
% base of the response.
%
% Input Variables
% Variable Name    Description
% pulse            The input pulse.
% threshold        The threshold for calculating the BW, normalized to
%                  [0,1]. For example, if set to 0.01, everything below
%                  0.01 in the response will be considered "not excited"
%                  or "not refocused".
% pulseType        'excitation' or 'refocusing'. If set to 'excitation',
%                  the BW of Mxy will be used starting from thermal 
%                  equilibrium. If set to 'refocusing', the BW of Mxy
% profileType      'mz' or 'mxy'. 

plotMin = -20;
plotMax = 20;
numPoints = 4000; % with a BW of 40 kHz, this gives 10 Hz resolution
vv = linspace(plotMin, plotMax, numPoints);
switch lower(pulseType)
    case {'excitation', 'inversion'}
        initMag = [0; 0; 1];
    case 'refocusing'
        initMag = [1; 0; 0];
    otherwise
        fprintf('Problem in CalcPulseBW: pulseType must be either excitation or refocusing. Currently set to %s. Aborting. \n', pulseType);
end

phaseAdjust = 1; % Apply a refocusing linear phase
response = PlotPulseFreqResponse(pulse, initMag, plotMin, plotMax, numPoints, profileType, [], 0, 1, phaseAdjust, 0);
response = response{1};

switch lower(profileType)
    case 'mz'
        switch lower(pulseType)
            case {'excitation', 'inversion'}
                baseline = 1;
            case 'refocusing'
                baseline = 0;
        end
    case 'mxy'
        switch lower(pulseType)
            case 'excitation'
                baseline = 0;
            case 'refocusing'
                baseline = 1;
        end
    otherwise
        fprintf('Problem in CalcPulseBW: profileType must be either mz or mxy. Currently set to %s. Aborting. \n', pulseType);
end        

if baseline==0
    idxBase = find(response>threshold, 1, 'first');
    BWBase  = abs(vv(idxBase))*2;
    idxFWHM = find(response>max(response)/2, 1, 'first');
    BWFWHM  = abs(vv(idxFWHM))*2;
%     halfResponse = response(end/2:end);
%     idxTop  = find(abs(halfResponse-halfResponse(1))>threshold, 1, 'first');
%     BWTop  = abs(vv(idxTop + numPoints/2))*2;
else
    idxBase = find(response<1-threshold, 1, 'first');
    BWBase  = abs(vv(idxBase))*2;
    idxFWHM = find(response<(1+min(response))/2, 1, 'first');
    BWFWHM  = abs(vv(idxFWHM))*2;
%     halfResponse = response(end/2:end);
%     idxTop = find(abs(halfResponse-halfResponse(1))>threshold, 1, 'first');
%     BWTop  = abs(vv(idxTop + numPoints/2))*2;
end