function fudgeFactor = CalcRefocusingFudgeFactor(pulse)
% Theoretically, when applying a 90-pulse for T ms in the presence of a 
% gradient G, a subsequent refocusing lobe of duration T/2 and gradient
% G (or any moment = G*T/2) are needed. However, in reality, this moment 
% may need to be slightly smaller or larger by a so-called "fudge factor", 
% close to unity, to achieve perfect refocusing. This function finds this 
% fudge factor (which may then be applied to either gradient or duration).


% T1 = 1e6;
% T2 = 1e6;
% freqMin = -
% CalcPulseFreqResponseWithRelax(pulse, T1, T2, freqMin, freqMax, numResPoints, [0; 0; 1]);