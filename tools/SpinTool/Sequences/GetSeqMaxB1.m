function maxB1 = GetSeqMaxB1(seq)
% Returns the maximal B1, in kHz, for a sequence object seq.
% For more information on sequences, consult ApplySequence.

maxB1 = 0;
numElements = numel(seq);
for idx=1:numElements
    curElement = seq{idx};
    if ispulse(curElement)
        maxB1 = max([maxB1 curElement.RFamp]);
    else
        cmd = curElement{1};
        switch lower(cmd)
            case 'hard'
                duration = 0.1;
                flipAngle = curElement{2};
                pulsePhase = curElement{3};
                pulse = PulseCreateHard(duration, flipAngle, pulsePhase);
                maxB1 = max([maxB1 pulse.RFamp]);
            case 'rect'
                duration = curElement{4};
                flipAngle = curElement{2};
                pulsePhase = curElement{3};
                pulse = PulseCreateHard(duration, flipAngle, pulsePhase);
                maxB1 = max([maxB1 pulse.RFamp]);
        end
    end
end