function db = kHz2DB(khz, tp90, power90db);
% The file outputs the pulse's maximal power in dB, given the power in kHz
% and the reference power of a 90 pulse.
%   khz         Input power, in khz.
%   tp90        90-pulse duration, in ms (e.g. 0.03 ms).
%   power90db   Power of the 90-pulse, in dB (e.g. 54 dB on the 500).

db = power90db + 20*log10(4*tp90*khz);