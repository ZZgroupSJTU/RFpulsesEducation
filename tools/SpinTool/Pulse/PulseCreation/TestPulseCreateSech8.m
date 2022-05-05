bw = 6;
b1 = 1.5;
T = 15;
N=512;

p = PulseCreateSech8(bw, 5.2, b1, T, N);
ph = PulseCreateSech(bw, 5.2, b1, T, N);
pc = PulseCreateChirp(-bw/2, bw/2, T, 1, N, 'z');
pc.RFamp = pc.RFamp./max(pc.RFamp)*b1;
%p2 = PulseShiftOffset(p, 10);
%p3 = PulseAddWithDelay(p, p2, 5);

PlotPulseFreqResponse({p, ph, pc}, [0; 0; 1], -5, 5, 300, 'mz');
%PlotPulseB1Performance(p, '1d', [0.2 1.5 40], 0, 'mz', [0; 0; 1]);