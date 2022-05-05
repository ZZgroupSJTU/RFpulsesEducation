function pulse = PulseScale(pulse,x)
% SYNTAX: pulse = PulseScale(pulse,x)
%
% Scales the duration of the pulse by an amount x, while appropriately
% scaling its amplitude by 1/x; 

if nargin<2
    x = 1;
end

if x<=0
    x = 1;
end

pulse.tp      = pulse.tp*x;
pulse.RFamp   = pulse.RFamp*1/x;
