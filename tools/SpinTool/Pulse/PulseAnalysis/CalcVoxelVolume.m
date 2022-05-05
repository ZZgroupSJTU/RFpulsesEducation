function [voxelVolume, signalVoxel] = CalcVoxelVolume(filenames, pulseDurations, flipAngles, VOISize, pathway, T1, TR, T2, TE)
% voxelVolume = CalcVoxelVolume(filename1, filename2, filename3, VOISize)
% 
% Given a cell array with three pulses (pulseArray), a voxel volume
% (VOISize, a 1x3 double vector) and which gradient axis.

plotRange = VOISize*2;

initMag = [0;0;1];
chemShift = 0;
B1Scaling = 1.0;
phaseAdjust = 0;
isPlot = 0;
numPoints = 141;

dV = plotRange(1)/numPoints * plotRange(2)/numPoints * plotRange(3)/numPoints; % in mm^3
dV = dV*0.001; % in cm^3

% ========================================================================
% Calculate the voxel volume of PRESS
% ========================================================================

%                            filename      pulseDuration      flipAngle      sliceThickness, gradientAxis)
p1 = PulseReadSiemensInclude(filenames{1}, pulseDurations(1), flipAngles(1), VOISize(1), 'z');
p2 = PulseReadSiemensInclude(filenames{2}, pulseDurations(2), flipAngles(2), VOISize(2), 'z');
p3 = PulseReadSiemensInclude(filenames{3}, pulseDurations(3), flipAngles(3), VOISize(3), 'z');

response1 = PlotPulseFreqResponse(p1,   initMag, -plotRange(1)/2, plotRange(1)/2, numPoints, 'flipangle', 'z', chemShift, B1Scaling, phaseAdjust, isPlot);
response2 = PlotPulseFreqResponse(p2,   initMag, -plotRange(2)/2, plotRange(2)/2, numPoints, 'flipangle', 'z', chemShift, B1Scaling, phaseAdjust, isPlot);
response3 = PlotPulseFreqResponse(p3,   initMag, -plotRange(3)/2, plotRange(3)/2, numPoints, 'flipangle', 'z', chemShift, B1Scaling, phaseAdjust, isPlot);

switch lower(pathway)
    case 'se'
        s1 = sind(response1{1}');
        s2 = sind(response2{1}'/2).^2;
        s3 = sind(response3{1}'/2).^2;

        s1Mat = repmat(s1, [1, numPoints, numPoints]);
        s2Mat = repmat(shiftdim(s2, -1), [numPoints, 1, numPoints]);
        s3Mat = repmat(shiftdim(s3, -2), [numPoints, numPoints, 1]);
        signalVoxel = s1Mat.*s2Mat.*s3Mat*exp(-TE/T2);
    case 'ste'
        s1 = sind(response1{1}');
        s2 = sind(response2{1}');
        s3 = sind(response3{1}');

        s1Mat = repmat(s1, [1, numPoints, numPoints]);
        s2Mat = repmat(shiftdim(s2, -1), [numPoints, 1, numPoints]);
        s3Mat = repmat(shiftdim(s3, -2), [numPoints, numPoints, 1]);

        signalVoxel = s1Mat.*s2Mat.*s3Mat/2*exp(-TE/T2/2)*exp(-TE/T1/2);
    case 'stress'
        s1 = sind(response1{1}');
        s2 = sind(response2{1}'/2).^2;
        s3 = sind(response3{1}'/2).^2;

        s1Mat = repmat(s1, [1, numPoints, numPoints]);
        s2Mat = repmat(shiftdim(s2, -1), [numPoints, 1, numPoints]);
        s3Mat = repmat(shiftdim(s3, -2), [numPoints, numPoints, 1]);

        s1B = sind(response1{1}');
        s2B = sind(response2{1}');
        s3B = sind(response3{1}');

        s1MatB = repmat(s1B, [1, numPoints, numPoints]);
        s2MatB = repmat(shiftdim(s2B, -1), [numPoints, 1, numPoints]);
        s3MatB = repmat(shiftdim(s3B, -2), [numPoints, numPoints, 1]);
        
        signalVoxel = s1Mat.*s2Mat.*s3Mat*exp(-TE/T2) + 0.5*s1MatB.*s2MatB.*s3MatB*exp(-TE/T2/2);
        
    otherwise
        error('Unrecognized pathway.');
end

voxelVolume = sum(signalVoxel(:)*dV);


