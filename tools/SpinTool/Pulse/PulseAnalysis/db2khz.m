function khz = db2khz(db, tp90, power90db);
% SYNTAX: khz = db2khz(db, tp90, power90db);
%
% The file outputs the pulse's maximal power in khz.
%   db          The input power, in dB, of the pulse.
%   tp90        90-pulse duration, in ms (e.g. 0.03 ms).
%   power90db   Power of the 90-pulse, in dB (e.g. 54 dB on the 500).


khz = 1/(4*tp90)*10^((db-power90db)/20);