function B=wurst(wurstParam, numSteps, fraction)
% SYNTAX: function B=wurst(wurstParam, numSteps);
% Returns a  "wurst ramp" function (of maximal size unity) of length N and
% rise exponential m. A typical value for m might be 24.
% fraction is an optional parameter between 0 and 1, indicating which
% fraction of the total number of steps the wurst function will occupy.

if nargin<3
    fraction = 1;
end

wurstSteps = round(numSteps*fraction);

p=[1:wurstSteps]-wurstSteps/2-0.5;
p=p./max(abs(p)).*pi./2;
B=1-abs(sin(p)).^wurstParam;

if fraction<1
    N = numSteps - wurstSteps;
    ceil(N/2)
    if N>0
        B = [zeros(1,ceil(N/2)), B];
    end
    if N>1
        B = [B, zeros(1, floor(N/2))];
    end
end