function pulse = OptimalControlPMExcitation(pulse, offsetVec, B1Vec, targetM, maxIter, epsilon, isOutput, maxB1)
% Designs an excitation pulse based on optimal control. Works only for
% excitation pulses (spins starting from +z). 
%
% Input Variables
% Variable Name     Units    Size   Description
% pulse             -        -      Initial guess pulse
% offsetVec         kHz      1xN    Vector of offsets
% B1Vec             scale    1xM    B1 scaling factor
% targetM           -        3xNxM  Target magnetization profile
% maxIter           -        #      Max. number of iterations
% epsilon           -        #      Gradient descent parameter
% isOutput          -        Bool   Set to 1 for feedback from function
% maxB1             kHz      #      Max. allowable amplitude of pulse

converged = 0;
idx = 0;
numFreqPts = size(targetM, 2);
numSteps = numel(pulse.RFamp);
freqMin = min(offsetVec);
freqMax = max(offsetVec);

M = zeros(numFreqPts, numSteps+1, 3);
% Lambda at the final time point = desired magnetization
L = zeros(numFreqPts, numSteps+1, 3);


% Keep track of the max. cost and the step at which it was encountered
maxIdx = 0;
maxCost = 0;

while (~converged)
    idx=idx+1;
    idxB1 = 1;
    for idxOffset=1:numFreqPts
        % Evolve initial states forward in time
        [Mx, My, Mz] = ApplyPulseDiagnostics(offsetVec(idxOffset), [0; 0; 0], [0; 0; 1], pulse);
        % Save M(t) for step IV
        M(idxOffset, :, 1) = Mx;
        M(idxOffset, :, 2) = My;
        M(idxOffset, :, 3) = Mz;
        % Evolve final target states backward in time
        finalState = squeeze(targetM(:, idxOffset, idxB1));
        [Lx, Ly, Lz] = ApplyPulseDiagnosticsInverse(offsetVec(idxOffset), [0; 0; 0], finalState, pulse);
        L(idxOffset, :, 1) = Lx;
        L(idxOffset, :, 2) = Ly;
        L(idxOffset, :, 3) = Lz;
    end
    % Calculate the field correction (AM only!)
    fieldCorrection = cross(M,L,3); % Cross along the 3rd dimension
    RF = pulse.RFamp.*exp(1i*pulse.RFphase);
    curCorrectionX = squeeze(fieldCorrection(:,:,1)); % Take the x-component of the cross product for current time
    curCorrectionX = sum(curCorrectionX,1); % Sum over offsets
    curCorrectionX = curCorrectionX(1:end-1);
    curCorrectionY = squeeze(fieldCorrection(:,:,2)); % Take the y-component of the cross product for current time
    curCorrectionY = sum(curCorrectionY,1); % Sum over offsets
    curCorrectionY = curCorrectionY(1:end-1);
    RF = RF - epsilon*(curCorrectionX + 1i*curCorrectionY);
    pulse.RFamp = abs(RF);
    pulse.RFphase = phase(RF);
    % Clip field
    pulse.RFamp(pulse.RFamp>maxB1) = maxB1;
    pulse.RFamp(pulse.RFamp<-maxB1) = -maxB1;
    % Check convergence
    curCost = Cost(pulse, targetM, freqMin, freqMax, B1Vec);
    if (isOutput), fprintf('Step %d cost %.7f (max %.7f at %d)  \n', idx, curCost, maxCost, maxIdx); end
    if (curCost>maxCost)
        maxCost = curCost;
        maxIdx = idx;
    end
    if (idx>maxIter)
        converged = 1;
    end
end

if (isOutput), fprintf('Max. cost was %.3f and was encountered at iteration %d\n', maxCost, maxIdx); end



% ========================================================================
%                            Auxiliary functions
% ========================================================================

function C = Cost(pulse, targetM, freqMin, freqMax, B1Vec)
% targetM: 3xN matrix, where targetM(:,i) is the target magnetization at
% freqAxis(i). Cost varies between -1 (worse) and +1 (best)
numFreqPts = size(targetM, 2);
numB1s = numel(B1Vec);
C = 0;
for idx=1:numB1s
    pulseTemp = pulse;
    pulseTemp.RFamp = pulseTemp.RFamp*B1Vec(idx);
    [Mx, My, Mz, ~, ~, ~] = CalcPulseFreqResponseWithRelax(pulseTemp, 1e6, 1e6, freqMin, freqMax, numFreqPts, [0; 0; 1], 0, 0);
    actualM = [Mx; My; Mz];
    C = sum(sum(actualM.*squeeze(targetM(:,:,idx))))/numFreqPts/numB1s;
end



