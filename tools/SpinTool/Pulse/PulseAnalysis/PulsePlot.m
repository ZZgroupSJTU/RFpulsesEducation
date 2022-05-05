function PulsePlot(pulse)
% The current function plots the pulse parameters of the input pulse object
% pulse can also be a cell array of pulses.

if iscell(pulse)
    numPulses = numel(pulse);
else
    numPulses = 1;
end

% Support up to 7 plots
colorOrder = [...
         0         0    1.0000;
         0    0.5000         0;
    1.0000         0         0;
         0    0.7500    0.7500;
    0.7500         0    0.7500;
    0.7500    0.7500         0;
    0.2500    0.2500    0.2500];

figure
ax(1)=subplot(2,4,1);
hold on;
ax(2)=subplot(2,4,2);
hold on;
ax(3)=subplot(2,4,3);
hold on;
ax(4)=subplot(2,4,4);
hold on;
ax(5)=subplot(2,4,5);
hold on;
ax(6)=subplot(2,4,6);
hold on;
ax(7)=subplot(2,4,7);
hold on;
ax(8)=subplot(2,4,8);
hold on;



for idx=1:numPulses
    if iscell(pulse)
        curPulse = pulse{idx};
        curColor = colorOrder(idx, :);
    else
        curPulse = pulse;
        curColor = colorOrder(1, :);
    end
    
    Nt = length(curPulse.RFamp);
    tt = linspace(0,curPulse.tp,Nt);

    axes(ax(1));
    plot(tt, curPulse.RFamp, 'color', curColor);
    curSAR(idx) = CalcSAR(curPulse, 'ref');
    title(['RF Amplitude. SAR (rel. to 1 ms 180): ', num2str(curSAR, 2)]);
    xlabel('ms');
    ylabel('kHz');

    axes(ax(2));
    plot(tt, curPulse.RFphase/2/pi, 'color', curColor);
    title('RF Phase');
    xlabel('ms');
    ylabel('2-\pi-rads.');

    axes(ax(3));
    Ot = (curPulse.RFphase([2:Nt]) - curPulse.RFphase([1:Nt-1]))/(curPulse.tp/Nt)/2/pi; 
    plot(linspace(0,curPulse.tp,Nt-1),Ot, 'color', curColor);
    xlabel('ms');
    title('Ot');
    ylabel('kHz');

    axes(ax(4));
    plot(tt, curPulse.Gx, 'color', curColor);
    title('Gx');
    xlabel('ms');
    ylabel('kHz/mm');

    axes(ax(5));
    plot(tt, curPulse.Gy, 'color', curColor);
    title('Gy');
    xlabel('ms');
    ylabel('kHz/mm');

    axes(ax(6));
    plot(tt, curPulse.Gz, 'color', curColor);
    title('Gz');
    xlabel('ms');
    ylabel('kHz/mm');

    axes(ax(7));
    plot(tt, real(curPulse.RFamp.*exp(1i*curPulse.RFphase)), 'color', curColor);
    title('Re.');
    xlabel('ms');
    ylabel('kHz');

    axes(ax(8));
    plot(tt, imag(curPulse.RFamp.*exp(1i*curPulse.RFphase)), 'color', curColor);
    title('Im.');
    xlabel('ms');
    ylabel('kHz');
end

linkaxes(ax, 'x');