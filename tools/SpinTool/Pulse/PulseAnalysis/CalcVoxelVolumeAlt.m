function [voxelVolumeSE, voxelVolumeSTE, voxelProfileSE, voxelProfileSTE, dV] = CalcVoxelVolumeAlt(filenames, pulseDurations, flipAngles, VOISize, plotRange, numSpins, B1Scaling)
% [voxelVolumeSE, voxelVolumeSTE, voxelProfileSE, voxelProfileSTE] 
%     = CalcVoxelVolumeAlt(filenames, pulseDurations, flipAngles, VOISize, B1Scaling)
% 
% Calculate the volume and profile of the voxel using the spin echo (SE) 
% and stimulated echo (STE) pathways. 

initMag = [0;0;1];
chemShift = 0;
if nargin<7
    B1Scaling = 1.0;
end
phaseAdjust = 0;
isPlot = 0;

dV = prod(plotRange)/prod(numSpins)*0.001; % in cm^3

numFlipAngles = size(flipAngles, 1);
numB1Values = numel(B1Scaling);

voxelVolumeSTE = zeros(numFlipAngles, numB1Values);
voxelVolumeSE = zeros(numFlipAngles, numB1Values);

% ========================================================================
% Calculate the voxel volume of PRESS
% ========================================================================

%                            filename      pulseDuration      flipAngles   sliceThickness, gradientAxis)
p1 = PulseReadSiemensInclude(filenames{1}, pulseDurations(1), 90.0,       VOISize(1), 'z');
p2 = PulseReadSiemensInclude(filenames{2}, pulseDurations(2), 90.0,       VOISize(2), 'z');
p3 = PulseReadSiemensInclude(filenames{3}, pulseDurations(3), 90.0,       VOISize(3), 'z');

p1Scaled = p1;
p2Scaled = p2;
p3Scaled = p3;

for idxB1=1:numB1Values
    for idxFlipAngle=1:numFlipAngles
        % Scale pulses amplitudes according to flip angles
        p1Scaled.RFamp = p1.RFamp/90.0*flipAngles(idxFlipAngle, 1);
        p2Scaled.RFamp = p2.RFamp/90.0*flipAngles(idxFlipAngle, 2);
        p3Scaled.RFamp = p3.RFamp/90.0*flipAngles(idxFlipAngle, 3);

        response1 = PlotPulseFreqResponse(p1Scaled, initMag, -plotRange(1)/2, plotRange(1)/2, numSpins(1), 'flipangle', 'z', chemShift, B1Scaling(idxB1), phaseAdjust, isPlot);
        response2 = PlotPulseFreqResponse(p2Scaled, initMag, -plotRange(2)/2, plotRange(2)/2, numSpins(2), 'flipangle', 'z', chemShift, B1Scaling(idxB1), phaseAdjust, isPlot);
        response3 = PlotPulseFreqResponse(p3Scaled, initMag, -plotRange(3)/2, plotRange(3)/2, numSpins(3), 'flipangle', 'z', chemShift, B1Scaling(idxB1), phaseAdjust, isPlot);

        % Calculate Spin Echo (SE)
        s1 = sind(response1{1}');
        s2 = sind(response2{1}'/2).^2;
        s3 = sind(response3{1}'/2).^2;

        s1Mat = repmat(s1,               [1, numSpins(2), numSpins(3)]);
        s2Mat = repmat(shiftdim(s2, -1), [numSpins(1), 1, numSpins(3)]);
        s3Mat = repmat(shiftdim(s3, -2), [numSpins(1), numSpins(2), 1]);

        voxelProfileSE{idxFlipAngle, idxB1} = s1Mat.*s2Mat.*s3Mat;
        voxelVolumeSE(idxFlipAngle, idxB1) = sum(voxelProfileSE{idxFlipAngle, idxB1}(:)*dV);

        % Calculate Spin Echo (STEVec)
        s1 = sind(response1{1}');
        s2 = sind(response2{1}');
        s3 = sind(response3{1}');

        s1Mat = repmat(s1,               [1, numSpins(2), numSpins(3)]);
        s2Mat = repmat(shiftdim(s2, -1), [numSpins(1), 1, numSpins(3)]);
        s3Mat = repmat(shiftdim(s3, -2), [numSpins(1), numSpins(2), 1]);

        voxelProfileSTE{idxFlipAngle, idxB1} = s1Mat.*s2Mat.*s3Mat/2;
        voxelVolumeSTE(idxFlipAngle, idxB1) = sum(voxelProfileSTE{idxFlipAngle, idxB1}(:)*dV);
    end
end