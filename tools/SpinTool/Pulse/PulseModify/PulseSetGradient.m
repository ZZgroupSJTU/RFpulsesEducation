function pulse = PulseSetGradient(pulse, gradient, gradientAxis)
% SYNTAX: 
%
%     pulse = PulseSetGradient(pulse, gradient, gradientAxis)
%
% Sets the gradient of the pulse object to a certain function along a given
% axis. The function (specified in the 'gradient' input) is interpolated
% if needed to fit the time points in the pulse object.
%
% Input Parameters
% Name           Units    Description
% pulse          -        Input pulse structure.
% gradient       kHz/mm   An 1xN vector of gradient values as a function
%                         of time (N=1 is possible). If N is different from
%                         the number of points in the pulse, the gradient
%                         shape is interpolated.
% gradientAxis   -        Set to either 'x', 'y', or 'z' (a character).
%
% Example: set the z-gradient to a constant, 2 kHz/mm
% pulse = PulseSetGradient(pulse, 2, 'z');

numPulseSteps = numel(pulse.RFamp);
numGradientSteps = numel(gradient);

if (numel(gradient)==1)
    gradientForm = ones(1,numPulseSteps).*gradient;
else
    if numel(gradient)==numel(pulse.RFamp)
        gradientForm = gradient;
    else
        % Interpolate gradient shape
        dt = pulse.tp/numPulseSteps;
        dtGrad = pulse.tp/numGradientSteps;
        pulseTimeAxis = [0:dt:(numPulseSteps-1)*dt];
        gradientTimeAxis = [0:dtGrad:(numGradientSteps-1)*dtGrad];
        gradientForm = interp1(gradientTimeAxis, gradient, pulseTimeAxis);
    end
end


switch (lower(gradientAxis))
    case 'x'
        pulse.Gx = gradientForm;
    case 'y'
        pulse.Gy = gradientForm;
    case 'z'
        pulse.Gz = gradientForm;
    otherwise
        disp('Error in SetGradient: gradient axis is not x, y or z! Aborting.');
end
