%% Chirp test
gammaHz=4275;
dimx.FOV=4;                % [cm]
dimx.FOVcenter=0;
dimx.FOVtype='z';
dimx.NFOV=256;
dimx.M0=1;
dimx.T2 = 10;
dimx.T1 = 10;

dimf.FOV=0;                % [cm]
dimf.FOVcenter=0;
dimf.FOVtype='cs';
dimf.NFOV=1;
dimf.M0=1;
dimf.T1 = 10;
dimf.T2 = 10;

density=0;
[M0,T1,T2,rcscell] = spin_system(density,dimx,dimf);


%%
pulse = PulseCreateChirpseg(-2.9, 2.9, 11.008, 8000,40,'z');
block{1}.tof=0;
block{1}.t=pulse.t;
block{1}.B1=pulse.RFamp.*pulse.RFphase/1000;
block{1}.Gx=pulse.Gx;
block{1}.Gy=pulse.Gy;
block{1}.Gz=pulse.Gz*10*1000/gammaHz;
    

[M,Mxyt,Mzt,t]=ConsoleSimulator(M0,M0,T1,T2,rcscell,block{1});


for k=1:1:5000;
N0=1;
Nend=256;
Nz=Nend-N0+1;
subplot(121);
quiver3(zeros(1,Nz),zeros(1,Nz),linspace(1,Nz,Nz),real(Mxyt(N0:Nend,k).'),imag(Mxyt(N0:Nend,k).'),zeros(1,Nz),2,'linewidth',1,'Color','b');

% subplot(121);PlotSpinsInCylinder(phase(Mxyt(N0:Nend,k)).',linspace(1,Nz,Nz),'r');

view([-40,14]);
% title(num2str(k));
axis off;

% subplot(122);quiver3(zeros(1,Nz),zeros(1,Nz),linspace(1,Nz,Nz),real(Mxyt(128+(N0:Nend),k).'),imag(Mxyt(128+(N0:Nend),k).'),zeros(1,Nz),2,'linewidth',1,'Color','b');
% subplot(122);PlotSpinsInCylinder(phase(Mxyt(128+(N0:Nend),k)).',linspace(1,Nz,Nz),'r');

% axis off;
title(num2str(k));
% view([-40,14]);
pause(0.5);
end