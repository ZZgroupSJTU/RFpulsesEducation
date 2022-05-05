%% 
% This is a RF pulse simulation demo using a Bloch simulator
% Created for ISMRM enducational talk of RF pulse in May 2022
% Zhiyong Zhang @ zhiyong.zhang@sjtu.edu.cn
% This is a demo code to show the Adiabatic pulses 
%% Part I: effective field and adiabatic condition
% Defined a 1D spin system
clc;clear;close all;
f=[-6000 0 6000];            % three spins with different frequencies
spins0.M0=ones(length(f),1);
spins0=Init_Spinsystem({'cs',f},spins0);
% create a chirp pulse structure using Assaf's SpinTool
BW=20;   % kHz
Tp=4;   % ms 
Np=10000; % Number of point
pulse=PulseCreateChirp(-BW/2,BW/2,Tp,Np);
pulse_Adabatic=pulse;
pulse_Adabatic.RFamp=BW/2*sqrt((1-linspace(-1,1,Np).^2));  %% To create a constant-amplitute effective field when frequency is sweeping
% create a sequence block with RF without a slice-selection gradient
block{1}=Generate_seqblock(0,...                                   %Reference frequency
                           Tp*1e-3,...                      %block's duration
                           length(pulse_Adabatic.RFamp),...               %block's number of point for simulation
                           pulse_Adabatic.RFamp.*exp(1i*pulse_Adabatic.RFphase),...  % B1 field (complex value)
                           0,...                                   % Gx gradient (support constant value and a shape)
                           0,...                                   % Gy gradient (support constant value and a shape)
                           0);                                     % Gz gradient (support constant value and a shape)
spins0.M=repmat([0 0 1],[length(f),1]);
[~,Mevolution1]=ConsoleSimulator(spins0,block{1});
%% show the effective field in B1 rotation frame
params.BGColor='k';
params.viewSeqblock=false;
params.SaveFormat='none';
params.viewBfield=3;  % 0->no Bfield only M, 1->Beff, 2->B0, B1, 3->B0 B1 Beff 4->Beff, M
% params.filename='Adabatic_B1B0Beff_pulse_demo';
params.Ntjump=20;
spin_idx=2;
params.B0t=Mevolution1.B0t(spin_idx,:);
figure();showMagnetizationWithSeqblock_Dynamic_B1RotateFrame(Mevolution1.Mt(spin_idx,:,:),block{1},params)
%% show B0,B1 field
params.BGColor='k';
params.viewSeqblock=false;
params.SaveFormat='none';
params.viewBfield=2;  % 0->no Bfield only M, 1->Beff, 2->B0, B1, 3->B0 B1 Beff 4->Beff, M
% params.filename='Adabatic_B1B0Beff_pulse_demo';
params.Ntjump=20;
spin_idx=2;
params.B0t=Mevolution1.B0t(spin_idx,:);
figure();showMagnetizationWithSeqblock_Dynamic_B1RotateFrame(Mevolution1.Mt(spin_idx,:,:),block{1},params)
%% show the effective field and spin magnetization at 0 frequency
params.BGColor='k';
params.viewSeqblock=false;
params.SaveFormat='none';
params.viewBfield=4;  % 0->no Bfield, 1->Beff, 2->B0, B1, 3->B0 B1 Beff 4->Beff, M
% params.filename='Adabatic_B1B0Beff_pulse_demo';
params.Ntjump=20;
spin_idx=2;
params.B0t=Mevolution1.B0t(spin_idx,:);
figure();showMagnetizationWithSeqblock_Dynamic_B1RotateFrame(Mevolution1.Mt(spin_idx,:,:),block{1},params)
%
%% Part II:  a chirp-180 pulse demo with a amplitute evolope of (1-cos^40(t)),
block{1}=Generate_seqblock(0,...                                   %Reference frequency
                           Tp*1e-3,...                      %block's duration
                           length(pulse_Adabatic.RFamp),...               %block's number of point for simulation
                           3*pulse_Adabatic.RFamp.*exp(1i*pulse_Adabatic.RFphase),...  % B1 field (complex value)
                           0,...                                   % Gx gradient (support constant value and a shape)
                           0,...                                   % Gy gradient (support constant value and a shape)
                           0);                                     % Gz gradient (support constant value and a shape)
seqblock=combineblocks(block);
spins0.M=repmat([0 0 1],[length(f),1]);
[spins1,Mevolution1]=ConsoleSimulator(spins0,block{1});
% show the effective field and spin magnetizations at w0=0kHz and w0=6kHz
params.BGColor='k';
params.viewSeqblock=false;
params.SaveFormat='none';
params.viewBfield=4;  % 0->no Bfield, 1->Beff, 2->B0, B1, 3->B0 B1 Beff 4->Beff, M
% params.filename='Adabatic_B1B0Beff_pulse_demo';
params.Ntjump=20;
spin_idx=[2 3];
params.B0t=Mevolution1.B0t(spin_idx,:);
figure();showMagnetizationWithSeqblock_Dynamic_B1RotateFrame(Mevolution1.Mt(spin_idx,:,:),block{1},params)

%
%% Part III:  adiabatic pulse insensitive to B1 inhomogeneities (by increasing the power of the pulse)
clc;clear;close all;
f=linspace(-8000,8000,256);            
spins0.M0=ones(length(f),1);
spins0=Init_Spinsystem({'cs',f},spins0);
%
BW=10;   % kHz
Tp=4;   % ms 
Np=1000; % Number of point
pulse=PulseCreateChirp(-BW/2,BW/2,Tp,Np);

block{1}=Generate_seqblock(0,...                                   %Reference frequency
                           Tp*1e-3,...                      %block's duration
                           length(pulse.RFamp),...               %block's number of point for simulation
                           pulse.RFamp.*exp(1i*pulse.RFphase),...  % B1 field (complex value)
                           0,...                                   % Gx gradient (support constant value and a shape)
                           0,...                                   % Gy gradient (support constant value and a shape)
                           0);                                     % Gz gradient (support constant value and a shape)

%%              
figure('Name','RFshape_demo','NumberTitle','off',...
        'color','k','position',[200 100 800 600/2]);
    
    
mp4name=fullfile(pwd,['Adabatic_demo_B1maxcreasing','.mp4']);
writerObj=VideoWriter(mp4name,'MPEG-4');
writerObj.FrameRate=4;
open(writerObj);
figure('Name','RFshape_demo','NumberTitle','off',...
        'color','k','position',[200 100 1600/2 600/2]);

MzTA=zeros(256,200);
factor=linspace(0.2,10,100);

for m=1:100
    spins0.M=repmat([0 0 1],[length(f),1]);
    block{1}.B1amp=pulse.RFamp*factor(m);
    spins=ConsoleSimulator(spins0,block{1});
    MzTA(:,m)=spins.M(:,3);
    xdata=linspace(0,pulse.tp,length(pulse.RFamp));
    ydata=block{1}.B1amp;
    subplot(121);plot(xdata,ydata,'y','linewidth',2);
    xmax=max(abs(xdata(:)));
    xmax=round(xmax*100)/100;
%     ymax=max(abs(ydata));
    ymax=4;
    ymax=round(ymax*100)/100;
    axis([0,xmax,-0.2*ymax ymax*2]);
    set(gca, 'color', 'k'); 
    set(gca,'XTick',linspace(0,xmax,3),'FontSize',12,'FontWeight','bold','XColor','w');
    set(gca,'xticklabels',{num2str(0),[num2str(xmax/2),' ms'], [num2str(xmax),' ms']});
    set(gca,'YTick',linspace(0,ymax,3),'FontSize',12,'FontWeight','bold','YColor','w');
    set(gca,'yticklabels',{num2str(0),[num2str(ymax/2),' kHz'], [num2str(ymax),' kHz']});


    xdata=f;
    ydata=spins.M(:,3);
    subplot(122);plot(xdata,ydata,'y','linewidth',2);
    xmax=max(abs(xdata(:)));
    xmax=round(xmax*100)/100;
    ymax=max(abs(ydata));
    ymax=round(ymax*100)/100;
    axis([-xmax,xmax,-1.5*ymax ymax*1.5]);
    set(gca, 'color', 'k'); 
    set(gca,'XTick',linspace(0,xmax,3),'FontSize',12,'FontWeight','bold','XColor','w');
    set(gca,'xticklabels',{[num2str(-xmax),' Hz'],num2str(0), [num2str(xmax),' Hz']});
    set(gca,'YTick',linspace(-ymax,ymax,3),'FontSize',12,'FontWeight','bold','YColor','w');
    set(gca,'yticklabels',{'-1' num2str(0),'1'});

    frame=getframe(gcf);
    writeVideo(writerObj,frame);
end        

close(writerObj);
%% plot the angles following with the B1 max (kHz)
B1max=linspace(0.2,10,100)*max(pulse.RFamp(:));
flipangle=180*acos(MzTA(end/2+1,1:100)/1)/pi;
figure('Name','RFshape_demo','NumberTitle','off',...
        'color','k','position',[200 100 800 600/2]);

xdata=B1max(:);
ydata=flipangle(:);
plot(xdata,ydata,'y','linewidth',2);
xmax=max(abs(xdata(:)));
xmax=round(xmax*100)/100;
ymax=max(abs(ydata));

ymax=round(ymax);
axis([0,xmax,-0.2*ymax ymax*1.2]);
set(gca, 'color', 'k'); 
set(gca,'XTick',linspace(0,xmax,3),'FontSize',12,'FontWeight','bold','XColor','w');
set(gca,'xticklabels',{num2str(0),[num2str(xmax/2),' kHz'], [num2str(xmax),' kHz']});
set(gca,'YTick',linspace(0,ymax,3),'FontSize',12,'FontWeight','bold','YColor','w');
set(gca,'yticklabels',{num2str(0),[num2str(90),'\circ'], [num2str(180),'\circ']});
