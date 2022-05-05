function [freqAxis,matMz] = PlotPulseProfiles1D(pulses,...
                                                T1,...
                                                T2,...
                                                freqMin,...
                                                freqMax,...
                                                numResPoints,...
                                                initMagState,...
                                                linPhaseCorrect,...
                                                phaseWrap,...
                                                plotResults,...
                                                axesHandles)

% Description: simulates and plots (if figure axes are given)
% the 1D responses of the Hadamard pulses. Only Mz is computed
% and displayed.
%
% Inputs
%
% Var. Name        Units   Description
% pulses           -       Structure containing RF pulses.
% T1               sec     Longitudinal relaxation const.
% T2               sec     Transverse relaxation const.
% freqMin          kHz     Min. freq to plot from.
% freqMax          kHz     Max. freq to plot to.
% numResPoints     -       Resolution of plots.
% initMagState     -       Initial magnetization, [Mx; My; Mz]
% linPhaseCorrect  -       Correct for linear phase of sim.?
% phaseWrap        -       Wrap the phase in [0,2pi]?
% plotResults      -       Plot results?
% axesHandles      -       If supplied, the function will
%                          plot onto these axes handles. If
%                          not, a special figure will be
%                          generated.
%
% Outputs
%
% Var. Name        Units   Description
% matMz            -       matMz(i,:) contains Mz as a 
%                          function of frequency, as obtained
%                          by simulating the i-th pulse.
% freqAxis         kHz     Frequency axis for plotting.

numPulses = length(pulses);
matMz = zeros(numPulses, numResPoints); % Pre-allocate
for numCurPulse = 1:numPulses
    [~, ~, Mz, ~, ~, freqAxis] ... 
        = CalcPulseFreqResponseWithRelax( ...
        pulses{numCurPulse}, ...
        T1*1000, ...
        T2*1000, ...
        freqMin, ...
        freqMax, ...
        numResPoints, ...
        initMagState, ...
        linPhaseCorrect, ...
        phaseWrap);
    matMz(numCurPulse,:) = Mz;
end    
hadamardMatrix = hadamard(numPulses);
% Plot results
if (plotResults == 1)
    if (~exist('axesHandles'))
        figure
        set(gcf,'Name','Simulated 1D excitation profiles');
        for idx=1:numPulses
            axesHandles(idx)=subplot(numPulses, 1, idx);
        end
    end
    for numCurPulse = 1:numPulses
        axes(axesHandles(numCurPulse));
        plot(freqAxis, matMz(numCurPulse, :));
        axis([min(freqAxis) max(freqAxis) -1 1]);
        xlabel('kHz');
        ylabel('M_z (a.u.)');
        title(['Pulse number: ',num2str(hadamardMatrix(numCurPulse,:))]);
    end
end