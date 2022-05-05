function pathways = CalculatePathwaySpoilingJ(k)
% Syntax: pathways = CalculatePathwaySpoilingJ(k)
%
% In the presence of N pulses and N spoiling gradients (one following
% each pulse), many coherence pathways will be formed. This function
% returns a list of all such coherence pathways, given the spoiling
% moments inbetween the pulses. The input k is a list (Nx3 matrix) of
% spoiling moments between the pulses. For example, for three pulses
%
%      _               _                
%     | |             | |                
% RF _| |_____________| |______________
%
%         ___________     __________
% Gx ____|    k1x    |___|    k2x   |_
%
%         ___________     __________
% Gy ____|    k1y    |___|    k2y   |_
%
%         ___________     __________
% Gz ____|    k1z    |___|    k2z   |_
%
%
% This will be described by the following input structure:
%
% k = [k1x  k1y  k1z] 
%     [k2x  k2y  k2z]   
%
% And will generate the following coherence pathways at its conclusion:
%
% Pathway     Coherence        Gradient Moment*
% 1           (0,1,1)          k1+k2
% 2           (0,1,0)          k1
% 3           (0,1,-1)         k1-k2
% 4           (0,0,1)          k2
% 5           (0,0,0)          0
% 6           (0,0,-1)         -k2 
% 7           (0,-1,1)         -k1+k2
% 8           (0,-1,0)         -k1  
% 9           (0,-1,-1)        -k1-k2  
% * - Subtraction is vectorial
%
% If we choose k2=k1 (vectorially), then we are only left with pathways
% #5 (longitudinal, undetectable), #3 (excited, flipped) and #7 (excited,
% flipped.). If we limit ourselves to pathways defined by a final coherence
% of -1, we are only left with pathway #3 (quadrature detection).
%
% See also: PrintCoherencePathways, PrintCoherencePathwaysJ,
%           CalculatePathwaySpoiling

numPulses = size(k,1);

if nargin<2
    flipAngles = 90*ones(1, numPulses);
end

flipAngles = flipAngles/180*pi; % Convert to radians

%              State           Pulse effects   Spoiling moment    Relative weight
pathways{1} = {'Longitudinal', '',             [0 0 0],           1}; % Start from thermal equilibrium


for idxPulse=1:numPulses
    pathwayCounter = 0;
    numPathways = numel(pathways);
    for idxPathway=1:numPathways
        pulseSeries = pathways{idxPathway}{2};
        moment = pathways{idxPathway}{3};
        weighting = pathways{idxPathway}{4};
        switch lower(pathways{idxPathway}{1})
            case 'longitudinal'
                % A longitudinal coherence can either be ... 
                % (1) Unaffected by the pulse, remain longitudinal and be
                %     unaffected by the next spoiler, OR get inverted
                %     by remain longitudinal
                pathwayCounter = pathwayCounter + 1;
                newPathways{pathwayCounter} = {'Longitudinal', [pulseSeries,'N'], moment, weighting*(cos(flipAngles(idxPulse)/2)^2-sin(flipAngles(idxPulse)/2)^2)};
                if (moment==0)
                    % (2) Excited
                    pathwayCounter = pathwayCounter + 1;
                    newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'E'], k(idxPulse,:), weighting*sin(flipAngles(idxPulse))};
                else
                    % (2) Excited with same moment
                    pathwayCounter = pathwayCounter + 1;
                    newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'E'], moment + k(idxPulse,:), weighting*sin(flipAngles(idxPulse))};
                    % (3) Excited with reversed moment
                    pathwayCounter = pathwayCounter + 1;
                    newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'E'], -moment + k(idxPulse,:), weighting*sin(flipAngles(idxPulse))};
                end
            case 'transverse'
                % A transverse coherence can either be ... 
                % (1) Unaffected by the pulse, remain transverse and be
                %     affected by the next spoiler
                pathwayCounter = pathwayCounter + 1;
                newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'N'], moment + k(idxPulse, :), weighting*cos(flipAngles(idxPulse)/2)^2};
                % (2) Affected by the pulse, get flipped and be affected 
                %     by the next spoiler
                pathwayCounter = pathwayCounter + 1;
                newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'F'], -moment + k(idxPulse, :), weighting*sin(flipAngles(idxPulse)/2)^2};
                % (3) Get stored
                pathwayCounter = pathwayCounter + 1;
                newPathways{pathwayCounter} = {'Longitudinal', [pulseSeries,'S'], moment, weighting*sin(flipAngles(idxPulse))};
        end
    end
    pathways = newPathways;
end

