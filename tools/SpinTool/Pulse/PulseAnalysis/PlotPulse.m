function figHandle = PlotPulse(pulse)
% The current function plots the pulse parameters of the input pulse object
% OLD FUNCTION, USE PULSEPLOT INSTEAD

Nt = length(pulse.RFamp);
tt = linspace(0,pulse.tp,Nt);

figHandle = figure;
ax(1) = subplot(2,4,1);
plot(tt, pulse.RFamp);
title('RF Amplitude');
xlabel('ms');
ylabel('kHz');

ax(2) = subplot(2,4,2);
plot(tt, pulse.RFphase/2/pi);
title('RF Phase');
xlabel('ms');
ylabel('2-\pi-rads.');

ax(3) = subplot(2,4,3);
Ot = (pulse.RFphase([2:Nt]) - pulse.RFphase([1:Nt-1]))/(pulse.tp/Nt)/2/pi; 
plot(linspace(0,pulse.tp,Nt-1),Ot);
xlabel('ms');
title('Ot');
ylabel('kHz');

ax(4) = subplot(2,4,4);
plot(tt, pulse.Gx*1000/42.576);
title('Gx');
xlabel('ms');
ylabel('mT/m');

ax(5) = subplot(2,4,5);
plot(tt, pulse.Gy*1000/42.576);
title('Gy');
xlabel('ms');
ylabel('mT/m');

ax(6) = subplot(2,4,6);
plot(tt, pulse.Gz*1000/42.576);
title('Gz');
xlabel('ms');
ylabel('mT/m');

ax(7) = subplot(2,4,7);
plot(tt, real(pulse.RFamp.*exp(1i*pulse.RFphase)));
title('Re.');
xlabel('ms');
ylabel('kHz');

ax(8) = subplot(2,4,8);
plot(tt, imag(pulse.RFamp.*exp(1i*pulse.RFphase)));
title('Im.');
xlabel('ms');
ylabel('kHz');

linkaxes(ax, 'x');