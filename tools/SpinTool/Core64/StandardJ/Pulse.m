classdef Pulse
    % Class for storing RF pulses
    
    properties
        dtAxis      % 1xN time axis (does not have to be uniform), in ms
        amp         % 1xN amplitude of RF pulse, in uT
        phase       % 1xN phase of RF pulse, in radians
        offset      % Central transmit frequency of pulse, in kHz
        grad        % 3xN matrix of gradients along three logical axes (PE, RO, SL) (mT/m)
    end
    
    properties (Dependent)
        B1          % 1xN complex vector describing B1(t) in the transverse plane
        timeAxis
        numSteps 
        dwellTime
        duration
    end
        
    
    methods
        function pulse = Pulse(varargin)
            % Constructor
            p = inputParser;
            p.addParameter('dwellTime', [], @(x) x>0); 
            p.addParameter('duration', [], @(x) x>0);
            p.addParameter('dtAxis', [], @(x) isvector(x') && prod(x>0) && isreal(x));
            p.addParameter('amp', 0, @(x) isvector(x') && prod(x>=0) && isreal(x));
            p.addParameter('phase', 0, @(x) isvector(x') && prod(x>=0) && isreal(x));
            p.addParameter('offset', 0, @(x) isreal(x));
            p.addParameter('gradX', [], @(x) isreal(x));
            p.addParameter('gradY', [], @(x) isreal(x));
            p.addParameter('gradZ', [], @(x) isreal(x));
            p.parse(varargin{:});
            pulse.amp = p.Results.amp;
            pulse.phase = p.Results.phase;
            pulse.offset = p.Results.offset;
            N = numel(pulse.amp);
            if ~isempty(p.Results.dtAxis)
                pulse.dtAxis = p.Results.dtAxis;
            else
                if ~isempty(p.Results.duration)
                    pulse.dtAxis = ones(1,pulse.numSteps)*p.Results.duration/pulse.numSteps; 
                else
                    error('wowowow');
                end
            end
            isEmptyGradWaveform = [isempty(p.Results.gradX), isempty(p.Results.gradY), isempty(p.Results.gradZ)];
            numGradSteps = [numel(p.Results.gradX), numel(p.Results.gradY), numel(p.Results.gradZ)];
            isNonEmptyGradients = sum(numGradSteps)>0;
            if sum(abs(diff(numGradSteps(numGradSteps~=0))))>0
                error('Gradient non-empty waveforms provided with different numbers of elements: [%d, %d, %d]', numGradSteps(1), numGradSteps(2), numGradSteps(3));
            end
            if isNonEmptyGradients
                if max(numGradSteps)~=N
                    error('Gradient non-empty waveforms have a different number of elements than RF amplitude waveform: [%d, %d, %d] vs. %d', numGradSteps(1), numGradSteps(2), numGradSteps(3), N);
                end
            else
                pulse.grad
            end
            if numGradSteps(1)==0, gradX = zeros(1, 
            numEmptyGradWaveforms = isempty(p.Results.gradX) + isempty(p.Results.gradY) + isempty(p.Results.gradZ);
            switch numel(numEmptyGradWaveforms)
                case 3
                    pulse.grad = zeros(3, N);
                case {0, 1}
                if numel(p.Results.gradX)==numel(p.Results.gradY) && numel(p.Results.gradX)==numel(p.Results.gradZ)
                    pulse.grad = (
                else
                    
                end
                    
                case 2
            end
            end
                
            
        end

        
        
        function pulse = SetDuration(pulse, duration)
            pulse.dtAxis = ones(1, pulse.numSteps)*(duration/pulse.numSteps);
        end
        
        function duration = get.duration(pulse)
            % Returns total pulse duration, in ms
            duration = sum(pulse.dtAxis);
        end
        
        function numSteps = get.numSteps(pulse)
            % Returns total number of steps in pulse
            numSteps = numel(pulse.amp);
        end
        
        function Plot(pulse)
            % Plots pulse behavior as a function of time
            figure;
            
            ax(1) = subplot(1,4,1);
            plot(pulse.timeAxis, pulse.amp);
            title('RF Amplitude');
            xlabel('ms');
            ylabel('uT');

            ax(2) = subplot(1,4,2);
            plot(pulse.timeAxis, pulse.RFphase);
            title('RF Phase');
            xlabel('ms');
            ylabel('rads.');

            ax(3) = subplot(1,4,3);
            Ot = diff(pulse.amp)./diff(pulse.timeAxis)/2/pi; 
            plot(pulse.timeAxis(1:end-1),Ot);
            xlabel('ms');
            title('Frequency');
            ylabel('kHz');

            ax(4) = subplot(1,4,4);
            plot(pulse.timeAxis, real(pulse.amp.*exp(1i*pulse.phase)));
            plot(pulse.timeAxis, imag(pulse.amp.*exp(1i*pulse.phase)), 'r');
            title('Real & Imaginary');
            xlabel('ms');
            ylabel('kHz');
            
            linkaxes(ax, 'x');            
        end
        
        function pulse = Copy(pulse)
        end
        
        function PlotFreqResponse(pulse)
        end
        
        
    end
    
end