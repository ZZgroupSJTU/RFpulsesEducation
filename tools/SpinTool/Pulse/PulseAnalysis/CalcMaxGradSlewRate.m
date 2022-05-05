function maxGradSlewRate = CalcMaxGradSlewRate(grad, duration)
% maxGradSlewRate = CalcMaxGradSlewRate(grad, duration)
%
% Returns maximal slew rate in mT/m/ms. Grad is either a vector in mT/m,
% or a pulse object (with gradients in kHz/mm). The duration is supplied
% in ms. 

if isstruct(grad)
    numSteps = numel(grad.RFamp);
    dt = grad.tp/numSteps; % ms
    Gx = grad.Gx;
    Gy = grad.Gy;
    Gz = grad.Gz;
    slewRateX = max(abs(diff(Gx)/dt));  % kHz/mm/ms
    slewRateY = max(abs(diff(Gy)/dt));  % kHz/mm/ms
    slewRateZ = max(abs(diff(Gz)/dt));  % kHz/mm/ms
    maxGradSlewRate = [slewRateX, slewRateY, slewRateZ]; % kHz/mm/ms
    maxGradSlewRate = maxGradSlewRate*1000/42.57; % mT/m/ms
else
    numSteps = numel(grad);
    dt = duration/numSteps; % ms
    grad = [0, grad, 0];
    slewRateReal = max(abs(diff(real(grad))/dt));  % mT/m/ms
    slewRateImag = max(abs(diff(imag(grad))/dt));  % mT/m/ms
    maxGradSlewRate = max(slewRateReal, slewRateImag);
    
end
