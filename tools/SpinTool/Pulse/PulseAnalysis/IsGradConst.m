function isConst = IsGradConst(pulse, axis)
% Checks if the gradient waveforms of a given pulse are constant in time.
% axis is an optional parameter, set to either 'x', 'y' or 'z' to indicate
% checking along a particular axis. If omitted, all axes will be checked.

isConst = true;

numStepsX = numel(pulse.Gx);
numStepsY = numel(pulse.Gy);
numStepsZ = numel(pulse.Gz);

gradFirst = [pulse.Gx(1) pulse.Gy(1) pulse.Gz(1)];
epsilon = 0.000001;

if nargin<2
    isCheckAxis = [1 1 1];
else
    switch lower(axis)
        case 'x'
            isCheckAxis = [1 0 0];
        case 'y'
            isCheckAxis = [0 1 0];
        case 'z'
            isCheckAxis = [0 0 1];
        otherwise
            error('Axis should be x, y, or z.');
    end
end

if isCheckAxis(1)
    numDifferentIndices = sum(abs(pulse.Gx - pulse.Gx(1))>epsilon);
    if numDifferentIndices>0
        isConst = false;
        return
    end
end

if isCheckAxis(2)
    numDifferentIndices = sum(abs(pulse.Gy - pulse.Gy(1))>epsilon);
    if numDifferentIndices>0
        isConst = false;
        return
    end
end

if isCheckAxis(3)
    numDifferentIndices = sum(abs(pulse.Gz - pulse.Gz(1))>epsilon);
    if numDifferentIndices>0
        isConst = false;
        return
    end
end
