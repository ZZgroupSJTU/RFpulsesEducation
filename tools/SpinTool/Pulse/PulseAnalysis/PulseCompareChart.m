function pulseData = PulseCompareChart(pulses, pulseNames, filename)
% Compares a cell array of pulses. Exports results to filename in Excel
% format.

if ~iscell(pulses)
    error('pulses should be a cell array of pulses.');
end

numPulses = numel(pulses);

if nargin<2
    pulseNames = {};
end

for idx=1:numPulses
    if ~isempty(pulseNames)
        pulseData.name{idx} = pulseNames{idx};
    else
        pulseData.name{idx} = 'No Name';
    end
    pulseData.duration(idx) = pulses{idx}.tp;
    pulseData.SARRef180(idx) = CalcSAR(pulses{idx}, 'ref');
    pulseData.SARRefMax(idx) = CalcSAR(pulses{idx}, 'rect');
    pulseData.maxB1(idx) = max(abs(pulses{idx}.RFamp));
    [~, FWHM] = CalcPulseBW(pulses{idx}, 0.1, 'inversion', 'mz');
    pulseData.FWHM(idx) = FWHM;
end

if nargin>2
    fid = fopen(filename, 'w');
    fprintf(fid, 'Filename \t Duration (ms) \t SAR (rel. 90) \t SAR (rel. rect) \t Peak B1 (kHz) \t BW (FWHM) \t BW per SAR \n');
    for idx=1:numPulses
        fprintf(fid, '%s \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f', ...
            pulseData.name{idx}, ...
            pulseData.duration(idx), ...
            pulseData.SARRef180(idx), ...
            pulseData.SARRefMax(idx), ...
            pulseData.maxB1(idx), ...
            pulseData.FWHM(idx), ...
            pulseData.FWHM(idx)/pulseData.SARRef180(idx));
        fprintf(fid, '\n');
    end
    fclose(fid);
end

