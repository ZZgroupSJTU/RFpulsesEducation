function pulse = PulseAddGradRamp(pulse, gradientAxis, gradRampTime)

dwellTime = pulse.tp/numel(pulse.RFamp);
numRiseSteps = ceil(gradRampTime/dwellTime);
truegradRampTime = numRiseSteps*dwellTime;
switch lower(gradientAxis)
    case 'x'
        G = pulse.Gx;
    case 'y'
        G = pulse.Gy;
    case 'z'
        G = pulse.Gz;
    otherwise
        error('gradientAxis = %s, should be x, y or z', gradientAxis);
end

dAmpRise = G(1)/numRiseSteps;
dAmpFall = G(end)/numRiseSteps;
pulse.tp = pulse.tp + 2*truegradRampTime;
pulse.RFamp = [zeros(1, numRiseSteps), pulse.RFamp, zeros(1,numRiseSteps)];
pulse.RFphase = [zeros(1, numRiseSteps), pulse.RFphase, zeros(1,numRiseSteps)];
switch lower(gradientAxis)
    case 'x'
        pulse.Gx = [0:dAmpRise:G(1)-dAmpRise, pulse.Gx, G(end)-dAmpFall:-dAmpFall:0];
        pulse.Gy = [zeros(1, numRiseSteps), pulse.Gy, zeros(1,numRiseSteps)];
        pulse.Gz = [zeros(1, numRiseSteps), pulse.Gz, zeros(1,numRiseSteps)];
    case 'y'
        pulse.Gx = [zeros(1, numRiseSteps), pulse.Gx, zeros(1,numRiseSteps)];
        pulse.Gy = [0:dAmpRise:G(1)-dAmpRise, pulse.Gy, G(end)-dAmpFall:-dAmpFall:0];
        pulse.Gz = [zeros(1, numRiseSteps), pulse.Gz, zeros(1,numRiseSteps)];
    case 'z'
        pulse.Gx = [zeros(1, numRiseSteps), pulse.Gx, zeros(1,numRiseSteps)];
        pulse.Gy = [zeros(1, numRiseSteps), pulse.Gy, zeros(1,numRiseSteps)];
        pulse.Gz = [0:dAmpRise:G(1)-dAmpRise, pulse.Gz, G(end)-dAmpFall:-dAmpFall:0];
end
