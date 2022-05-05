clc;
clear all;
close all;

%%

refpw=10e-6;
reftpwr=58;
gamaHz=4257.4;

%%

p1_rf=init_rf('zz_SLRtest');
p1_rf.duration=0.01;
p1_rf.pts=1000;
p1_rf.flip=90;
p1_rf.bandwidth=0.5;                          % in [kHz] : RF pulse bandwidth

SLRParams=init_SLRparam(p1_rf.duration,p1_rf.bandwidth,0.0001,0.0001,'equiripple','excitation',p1_rf.pts,p1_rf.flip);

[a, b] = calc_SLRpolyn(SLRParams);

[thetalong, philong] = inverseSLR(a,b);

p1_rf=calc_SLRrf(p1_rf,thetalong,philong,refpw,reftpwr);

saveRFfile(p1_rf,'/home/Zhiyongi/Program/NMRSimulator/shapelib');