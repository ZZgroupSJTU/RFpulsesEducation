SLR_RF.pulname='zz_SLR90';
SLR_RF.duration=0.004;
SLR_RF.pts=1000;
SLR_RF.flip=90;
SLR_RF.bandwith=1.667;
SLR_RF.inripple=0.001;
SLR_RF.outripple=0.001;
SLR_RF.filtertype='leastsquare';
SLR_RF.pulsetype='excitation';
SLR_RF.multiws=1;
SLR_RF.multiwp=1;

SLR_RF.multiband='false';
p1_rf=CreateSLR_RF(SLR_RF);
%%
figure();plot(p1_rf.amp)
% figure();plot(p1_rf.phase)