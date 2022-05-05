function PlotPulseFMTrajectory(pulse, offset)
% Plots, in the Frequency Modulated (FM) frame, the trajectory of the 
% given pulse.

Nt = numel(pulse.RFamp);
Ot = (pulse.RFphase([2:Nt]) - pulse.RFphase([1:Nt-1]))/(pulse.tp/Nt)/2/pi; 
Ot(end+1)=Ot(end);

wx = pulse.RFamp;
wz = offset + Ot;

figure
plot(wx,wz);
xlabel('kHz');
ylabel('kHz');