classdef SequenceElement
    % A PulseSequence element. 
    %   name        String of element's name. Can be left empty.
    %   type        String describing element type. Possible options (case
    %               insensitive): 'Pulse', 'Delay', 'CSIEncode', 'Purge',
    %               'PurgeMoment', 'Acquire', 'FreqEncode', 'HardPulse',
    %               'ConstPulse'
    %   logicalVol  VDIVolume describing the relation between the logical
    %               and physical gradient axes
    %   maxSlewRate 
    %   maxGrad
    %   maxB1
    %   
    %   Element-specific properties:
    %   
    %   Pulse     
    %     pulse     A Pulse object
    %
    %   Delay
    %     duration     (ms)
    %     numSteps    
    %     
    %   CSIEncode
    %     encodeType   Case-insensitive string. Determines how the CSI
    %                  gradients will be calculated: 
    %                  'minimumTime': Maximal gradient slew rates to 
    %                                 minimize the total CSI encoding time
    %                  'fixedTime':   Use the 
    %     FOV          1x3 vector of FOV of CSI encoding (mm)
    %     step         Current phase encoding step, from [0,0,0] to 
    %                  matrixSize-[1 1 1].
    %     matrixSize   1x3 vector of matrix size
    %     maxSlewRate  Maximal gradient slew rate, assumed equal along all
    %                  axes. Only used if encodeType is set to minimumTime
    %                  (mT/m/ms)
    %     duration     ms
    %     
    %  Acquire
    %     SW
    %     dwellTime
    %
    %  FreqEncode
    %     FOV
    %     BWPerPixel
    %     acqAxis
    %
    %  Purge
    %     grad        1x3 vector of gradient intensities (mT/m)
    %     duration
    %
    %  PurgeMoment
    %     moment      1x3 vector of gradient moments (m^(-1))
    %     duration
    %
    %  HardPulse
    %     numSteps
    %     duration
    %     flipAngle
    %     pulsePhase
    %
    %  ConstPulse
    %     numSteps
    %     duration
    %     pulseAmp
    %     pulsePhase
    
    properties
        name
        type
        prop
    end
    
    properties (Dependent)
        duration
    end
    
    methods
        function pulse = GetPulse(element)
            % GetPulse  Returns the Pulse associated with the sequence element
            switch lower(element.type)
                case 'pulse'
                    pulse = element.prop.pulse;
                case 'hardpulse'
                    pulse = 
                otherwise
                    error('Unrecognized element type %s\n', element.type);
            end
        end
    end
end