classdef PulseSequence
    % Class for storing pulse sequences.
    %   Sequence is stored as a structured array of 'elements', which
    %   consist of pulses, delays, CSI encoding gradients and so forth.
    %   All elements share several fields, and then each element type
    %   will have its own unique field as well.
    %   
    %   
    properties
        element
        maxGrad
        maxSlewRate
        maxB1
    end
    
    properties (Dependent)
        totalTime
    end
end