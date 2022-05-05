function pathways = CalculatePathwaySpoiling(k, flipAngles)
% Syntax: pathways = CalculatePathwaySpoiling(k, flipAngles)
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
%         ___________
% Gx ____|    k1x    |_________________
%
%                          __________
% Gy _____________________|    k2y   |_
%
%                                     
% Gz __________________________________
%
%
% This will be described by the following input structure:
%
% k = [k1x  0   0  ] 
%     [ 0  k2y  0  ]   
%
% And will generate the following coherence pathways at its conclusion:
%
% Pathway     State            Pulse 1?        Pulse 2?      Gradients*
% 1           Transverse       Excite          Flip          k2-k1 
% 2           Transverse       Excite          No effect     k2+k1
% 3           Longitudinal     Excite          Store         k1
% 4           Longitudinal     No effect       Excite        k2
% * - Subtruction is vectorial
%
% This will be returned as the following cell array list:
% pathways{1} = {'Transverse', 'EF', k2-k1, [-1 1],  0, [-1 2]}
% pathways{2} = {'Transverse', 'EN', k2+k1, [-1 -1], 1, [-1 0]} 
% etc ...
% The general format is:
% {Type of coherence, effect of pulses, gradient moment, p-coherence, flip angle weighting, p-coherence delta}
% The p-coherence is a vector having the same number of elements as the 
% number of RF pulses, indicating which coherence (+1,-1,0) the spins are
% at each stage throughout the sequence for this particular pathway.
% The p-coherence delta is a vector which shows the "deltas" of the 
% p-coherence vector. 
%
% The 'flipAngles' input is optional, and will be set to all 90s if
% not used. It contains the flip angles of the different pulses, for
% calculating the weighting of the pathway. Angles are in degrees.
%
% The rules for coherence pathway creation are detailed in Haacke.
%
% See also: PrintCoherencePathways

numPulses = size(k,1);

if nargin<2
    flipAngles = 90*ones(1, numPulses);
end

flipAngles = flipAngles/180*pi; % Convert to radians

%              State           Pulse effects   Spoiling moment    Spoiling Moment Matrix        Relative weight    RF phase matrix
pathways{1} = {'Longitudinal', '',             [0 0 0],           []                            1,                 []}; % Start from thermal equilibrium


for idxPulse=1:numPulses
    pathwayCounter = 0;
    numPathways = numel(pathways);
    for idxPathway=1:numPathways
        pulseSeries = pathways{idxPathway}{2};
        moment = pathways{idxPathway}{3};
        momentMatrix = pathways{idxPathway}{4}; 
        weighting = pathways{idxPathway}{5};
        RFMatrix = pathways{idxPathway}{6};
        switch lower(pathways{idxPathway}{1})
            case 'longitudinal'
                % A longitudinal coherence can either be ... 
                % (1) Unaffected by the pulse, remain longitudinal and be
                %     unaffected by the next spoiler, OR get inverted
                %     by remain longitudinal
                pathwayCounter = pathwayCounter + 1;
                newPathways{pathwayCounter} = {'Longitudinal', [pulseSeries,'N'], moment, [momentMatrix, 0], weighting*(cos(flipAngles(idxPulse)/2)^2-sin(flipAngles(idxPulse)/2)^2), [RFMatrix, 0]};
                if (moment==0)
                    % (2) Excited
                    pathwayCounter = pathwayCounter + 1;
                    newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'E'], k(idxPulse,:), [momentMatrix, 1], weighting*sin(flipAngles(idxPulse)), [RFMatrix, 1]};
                else
                    % (2) Excited with same moment
                    pathwayCounter = pathwayCounter + 1;
                    newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'E'], moment + k(idxPulse,:), [momentMatrix, 1], weighting*sin(flipAngles(idxPulse)), [RFMatrix, 1]};
                    % (3) Excited with reversed moment
                    pathwayCounter = pathwayCounter + 1;
                    newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'E'], -moment + k(idxPulse,:), [-momentMatrix, 1], weighting*sin(flipAngles(idxPulse)), [-RFMatrix, 1]};
                end
            case 'transverse'
                % A transverse coherence can either be ... 
                % (1) Unaffected by the pulse, remain transverse and be
                %     affected by the next spoiler
                pathwayCounter = pathwayCounter + 1;
                newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'N'], moment + k(idxPulse, :), [momentMatrix, 1], weighting*cos(flipAngles(idxPulse)/2)^2, [RFMatrix, 0]};
                % (2) Affected by the pulse, get flipped and be affected 
                %     by the next spoiler
                pathwayCounter = pathwayCounter + 1;
                newPathways{pathwayCounter} = {'Transverse', [pulseSeries,'F'], -moment + k(idxPulse, :), [-momentMatrix, 1], weighting*sin(flipAngles(idxPulse)/2)^2, [-RFMatrix, 2]};
                % (3) Get stored
                pathwayCounter = pathwayCounter + 1;
                newPathways{pathwayCounter} = {'Longitudinal', [pulseSeries,'S'], moment, [momentMatrix, 0], weighting*sin(flipAngles(idxPulse)), [RFMatrix -1]};
        end
    end
    pathways = newPathways;
end

