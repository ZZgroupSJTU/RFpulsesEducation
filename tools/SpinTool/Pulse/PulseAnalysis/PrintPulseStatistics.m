function PrintPulseStatistics(pulseArray)

if ~iscell(pulseArray)
    if ~isstruct(pulseArray)
        fprintf('Error in PrintPulseStatistics: input must be a pulse or a cell array of pulses. Aborting!\n');
        return
    end
    pulseArray = {pulseArray};
end

numPulses = numel(pulseArray);
pulseHard180 = PulseCreateHard(500, 180, 0, 50);
for idx=1:numPulses
    duration(idx) = pulseArray{idx}.tp;
    numSteps(idx) = numel(pulseArray{idx}.RFamp);
    maxB1(idx) = max(pulseArray{idx}.RFamp);
    SAR(idx) = CalcSAR(pulseArray{idx}, pulseHard180);
    gradient(idx) = max(pulseArray{idx}.Gz)*1000/GetGyromagneticRatio('1h');
end

fprintf('Pulse Statistics:\n');
fprintf('Pulse #:           '); fprintf('%d\t\t', [1:numPulses]); fprintf('\n');
fprintf('Duration (ms):     '); fprintf('%.2f\t\t', duration); fprintf('\n');
fprintf('Peak B1 (kHz):     '); fprintf('%.2f\t\t', maxB1); fprintf('\n');
fprintf('SAR*:              '); fprintf('%.2f\t\t', SAR); fprintf('\n');
fprintf('Num. Steps:        '); fprintf('%d\t\t', numSteps); fprintf('\n');
fprintf('z-Gradient (mT/m): '); fprintf('%.2f\t\t', gradient); fprintf('\n');
fprintf('* - Relative to SAR of a 0.5 ms 180 pulse\n');