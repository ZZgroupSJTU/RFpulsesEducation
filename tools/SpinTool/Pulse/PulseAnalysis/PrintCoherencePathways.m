function PrintCoherencePathways(x, thresholds, inclusion, voxelSize)
% SYNTAX: PrintCoherencePathways(x, thresholds, inclusion, voxelSize)
%
% To be used in conjunction with CalculatePathwaySpoiling's output. 
% Prints the difference coherence pathways; their moments; their state
% (longitudinal/transverse); the effects of the pulses; the flip angle
% weighting; and the spoiling of the signal for the given
% voxel size, assuming a homogeneous voxel (but not taking into
% account the flip angle weighting).
% 
% Input Variables
% Variable Name  Description
% x              x is the output of CalculatePathwaySpoiling.
% thresholds     Optional. A 2-element vector of [low, high] moments.
%                Anything outside this range will not be printed.
%                Defaults to [0 1e6] ("everything")
% inclusion      Optional. A string equal to 'longitudinal', 'transverse'
%                or anything else, and will show the listed component.
%                Defaults to 'transverse'.
% voxelSize      Optional. A 1x3 vector containing the VOI/voxel/ROI
%                size, in mm. Assumes spoiling moments are in m^(-1).
%                Defaults to [10 10 10] (typical voxel size).
%
% See also: CalculatePathwaySpoiling

numPathways = numel(x);
if nargin<2
    thresholds = [0 1e6];
end

if nargin<3
    inclusion = 'transverse';
end
if ~(strcmpi(inclusion,'transverse') || strcmpi(inclusion, 'longitudinal'))
    inclusion = 'transverse';
end

if nargin<4
    voxelSize = [10 10 10];
end

fprintf('Calculated Pathways: \n');
fprintf('Results are presented in table format: \n');
fprintf('   Pathway #:  Type (transverse/longitudinal), Pulse effects, k-space Moment, Attenuation, Signal Relative Amp.\n');
for idx=1:numPathways
    % x{k}{1} - Pulse type, string (longitudinal/transverse)
    % x{k}{2} - Pulse effects, string (e.g. 'EFF' - Excite flip flip)
    % x{k}{3} - Total gradient moment
    % x{k}{4} - Coherence level during the sequence for this particular pathway
    % x{k}{5} - Signal weighting (1: Full signal, no attenuation)
    % x{k}{6} - Difference of coherence levels during sequence for this particular pathway
    moment = x{idx}{3};
    momentInRange = ((norm(moment)<=thresholds(2)) && (norm(moment)>=thresholds(1)));
    if (strcmpi(x{idx}{1}, inclusion)) && momentInRange
        signalWeight = sinc(pi*moment(1)*voxelSize(1)*0.001)*sinc(pi*moment(2)*voxelSize(2)*0.001)*sinc(pi*moment(3)*voxelSize(3)*0.001);
        if (norm(moment)<1) % Note (with asterixes) moments that are very small
            fprintf('*** Pathway %d:  \t %s  \t  %s  [%.1f %.1f %.1f]    \t(%.3f) \t%.3f *** \n', idx, x{idx}{1}, x{idx}{2}, x{idx}{3}(1), x{idx}{3}(2), x{idx}{3}(3), x{idx}{5}, signalWeight);
        else
            fprintf('Pathway %d:  \t %s  \t  %s  [%.1f %.1f %.1f] \t(%.3f) \t%.3f \n', idx, x{idx}{1}, x{idx}{2}, x{idx}{3}(1), x{idx}{3}(2), x{idx}{3}(3), x{idx}{5}, signalWeight);
        end
    end
end
fprintf('\n');