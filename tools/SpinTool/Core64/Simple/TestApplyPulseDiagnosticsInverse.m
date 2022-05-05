clear all
close all
clc

isRecompile = 0;
if (isRecompile)
    mex ApplyPulseDiagnostics.c
    mex ApplyPulseDiagnosticsInverse.c
end

pulse = PulseCreateSinc(4, 4, 8, 90);

[Mx, My, Mz] = ApplyPulseDiagnostics(0, [0; 0; 0], [0; 0; 1], pulse);
M = [Mx(end); My(end); Mz(end)];
[MxInv, MyInv, MzInv] = ApplyPulseDiagnosticsInverse(0, [0; 0; 0], M, pulse);

figure
subplot(3,1,1)
plot(Mx);
hold
plot(MxInv, 'r');

subplot(3,1,2)
plot(My);
hold
plot(MyInv, 'r');

subplot(3,1,3)
plot(Mz);
hold
plot(MzInv, 'r');
