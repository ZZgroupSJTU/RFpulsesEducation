function R = CalcPulseSweepRate(pulse)
% Calculates the sweep rate of the pulse as a function of time. The 
% sweep rate is defined as the derivative of the pulse's instantaneous
% frequency, and is given in kHz^2.

Nt = length(pulse.RFamp);
dt = pulse.tp/Nt;

% Instantaneous frequency, in kHz (without the 2pi!)
Ot = diff(pulse.RFphase)/(2*pi*dt); 

R = diff(Ot)/dt;
R = [R(1) R R(end)];