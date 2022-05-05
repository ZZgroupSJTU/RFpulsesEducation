function [FWHM, TB] = PulseCalcInversionBW(pulse)

response = PlotPulseFreqResponse(pulse, [0;0;1], -20, 20, 1001, 'mz');
Mz = response{1};
minMz = min(Mz);
halfMax = (1+minMz)/2;


