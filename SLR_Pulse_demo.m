%% 
% This is a RF pulse simulation demo using a Bloch simulator
% Created for ISMRM enducational talk of RF pulse in May 2022
% Zhiyong Zhang @ zhiyong.zhang@sjtu.edu.cn
%% This is a demo code to show the SLR pulses     
clc; clear;close all;
%% Defined a 1D spin system
Nf=1000;                                   
f=linspace(-5000,5000,Nf);          %[Hz]
spins0.M0=ones(Nf,1);
spins0=Init_Spinsystem({'cs',f},spins0);
%% definite a RF profile information and create a pulse structure using Assaf's SpinTool
pulse_duration=0.004;   %[ms]
pulse_number_of_point=400;
pulse_number_of_lobe=6;
pulse_flipangle=180;
pulse_initphase=0;
pulse=PulseCreateSinc(pulse_number_of_lobe,...
                        pulse_duration*1e3,...
                        pulse_number_of_point,...
                        pulse_flipangle,...
                        pulse_initphase);

%% Create a SLR pulse

SLR_RF.pulname='zz_SLR90';
SLR_RF.duration=pulse_duration;
SLR_RF.pts=400;
SLR_RF.flip=pulse_flipangle;
SLR_RF.bandwith=6;
SLR_RF.inripple=0.00001;
SLR_RF.outripple=0.00001;
% SLRpulsetype='alpha' ,  'excitation' ,  'inversion',   'refocusing'
% SLRfiltertype= 'leastsquare' , 'equiripple' ,  'complexequiripple',  'quadraticequiripple' 
SLR_RF.filtertype='complexequiripple';
SLR_RF.pulsetype='inversion';
SLR_RF.multiws=1;
SLR_RF.multiwp=1;
SLR_RF.multiband='false';

p1_rf=CreateSLR_RF(SLR_RF);
pulseSLR.tp=p1_rf.duration*1e3;
pulseSLR.RFamp=p1_rf.amp;
pulseSLR.RFphase=p1_rf.phase;
pulseSLR.Gx=0;
pulseSLR.Gy=0;
pulseSLR.Gz=0;
%%

spins0.M=repmat([0 0 1],[Nf,1]);
block{1}=Generate_seqblock(0,...                             %Reference frequency
                           pulse.tp*1e-3,...                    %block's duration
                           length(pulse.RFamp),...             %block's number of point for simulation
                           pulse.RFamp.*exp(1i*pulse.RFphase),...  % B1 field (complex value)
                           0,...                                 % Gx gradient (support constant value and a shape)
                           0,...                                 % Gy gradient (support constant value and a shape)
                           0);                                   % Gz gradient (support constant value and a shape)
                       
spins1=ConsoleSimulator(spins0,block{1});

block{2}=Generate_seqblock(0,...                             %Reference frequency
                           pulseSLR.tp*1e-3,...                    %block's duration
                           length(pulseSLR.RFamp),...             %block's number of point for simulation
                           pulseSLR.RFamp.*exp(1i*pulseSLR.RFphase),...  % B1 field (complex value)
                           0,...                                 % Gx gradient (support constant value and a shape)
                           0,...                                 % Gy gradient (support constant value and a shape)
                           0);                                   % Gz gradient (support constant value and a shape)
spins2=ConsoleSimulator(spins0,block{2});

%% Show the 180 inversion comparison (SLR vs SINC)
figure('Name','RFshape_demo','NumberTitle','off',...
        'color','k','position',[200 100 1600/2 600/2]);

axes1=subplot(121);
t=block{1}.t;
plot(t,block{1}.B1amp,'r','linewidth',2); 
hold on;plot(block{2}.t,block{2}.B1amp,'y','linewidth',2);
axis([0,t(end),-0.2 4]);
set(axes1, 'color', 'k');  
set(axes1,'XTick',linspace(0,t(end),3),'FontSize',12,'FontWeight','bold','XColor','w');
set(axes1,'xticklabels',{num2str(0),[num2str(t(end)/2*1e3),' ms'], [num2str(t(end)*1e3),' ms']});
set(axes1,'YTick',[-0.5 0  0.5 1 1.5 2 2.5],'FontSize',12,'FontWeight','bold','YColor','w');
text(t(end)/3,3.6, 'SINC pulse \pi ', 'FontSize',10,'color','r')
text(t(end)/3,3.3, 'SLR pulse \pi ', 'FontSize',10,'color','y')
 
axes2=subplot(122);
% plot(f,abs(spins1.M(:,1)+1i*spins1.M(:,2)),'r','linewidth',2);
% hold on;plot(f,abs(spins2.M(:,1)+1i*spins2.M(:,2)),'y','linewidth',2);

% plot(f,spins1.M(:,2),'r','linewidth',2);
% hold on;plot(f,spins2.M(:,2),'y','linewidth',2);

plot(f,spins1.M(:,3),'r','linewidth',2);
hold on;plot(f,spins2.M(:,3),'y','linewidth',2);

f_bandwith=f(end)-f(1);

text(-f_bandwith/2.5,1.7, 'SINC pulse inversion profile (M_z)', 'FontSize',10,'color','r')
text(-f_bandwith/2.5,1.45, 'SLR pulse inversion profile (M_z)', 'FontSize',10,'color','y')
    
axis([-f_bandwith/2 f_bandwith/2 -1.2,2]);
%     title('Excited Magnetization(Mxy)','color','y');
set(axes2, 'color', 'k');  
set(axes2,'XTick',linspace(-f_bandwith/2,f_bandwith/2,3),'FontSize',12,'FontWeight','bold','XColor','w');
set(axes2,'xticklabels',{num2str(-round(f_bandwith/2*1e-2)/10),'0', [num2str(round(f_bandwith/2*1e-2)/10),' kHz']});
set(axes2,'FontSize',12,'FontWeight','bold','YColor','w');

%% Excitation comparison (SLR vs SINC) for increasing flip angle
clear spins0;

Nf=256;                                   
f=linspace(-2000,2000,Nf);          %[Hz]
spins0.M0=ones(Nf,1);
spins0=Init_Spinsystem({'cs',f},spins0);
% sinc pulse
pulse_duration=0.004;   %[ms]
pulse_number_of_point=400;
pulse_number_of_lobe=3;
pulse_flipangle=90 ;
pulse_initphase=0;
pulse=PulseCreateSinc(pulse_number_of_lobe,...
                        pulse_duration*1e3,...
                        pulse_number_of_point,...
                        pulse_flipangle,...
                        pulse_initphase);
% SLR using Assaf's SpinTool
pulseSLR = PulseCreateSLR90(400, 1, 10, 'linearR6', 0);   
pulseSLR.RFamp=pulseSLR.RFamp*pulseSLR.tp/pulse.tp; 
pulseSLR.tp=pulse.tp;

mp4name=fullfile(pwd,['SLR_demo_TAcreasing','.mp4']);
writerObj=VideoWriter(mp4name,'MPEG-4');
writerObj.FrameRate=4;
open(writerObj);
figure('Name','RFshape_demo','NumberTitle','off',...
        'color','k','position',[200 100 1600/2 600/2]);


for m=1:4:180 

    block{1}=Generate_seqblock(0,...                             %Reference frequency
                               pulse.tp*1e-3,...                    %block's duration
                               length(pulse.RFamp),...             %block's number of point for simulation
                               m/90*pulse.RFamp.*exp(1i*pulse.RFphase),...  % B1 field (complex value)
                               0,...                                 % Gx gradient (support constant value and a shape)
                               0,...                                 % Gy gradient (support constant value and a shape)
                               0);                                   % Gz gradient (support constant value and a shape)
    spins1=ConsoleSimulator(spins0,block{1});

    block{2}=Generate_seqblock(0,...                             %Reference frequency
                               pulseSLR.tp*1e-3,...                    %block's duration
                               length(pulseSLR.RFamp),...             %block's number of point for simulation
                               m/90*pulseSLR.RFamp.*exp(1i*pulseSLR.RFphase),...  % B1 field (complex value)
                               0,...                                 % Gx gradient (support constant value and a shape)
                               0,...                                 % Gy gradient (support constant value and a shape)
                               0);                                   % Gz gradient (support constant value and a shape)
    spins2=ConsoleSimulator(spins0,block{2});
   

    axes1=subplot(121);
    t=block{1}.t;
    plot(t,block{1}.B1amp*1.5,'r','linewidth',2); 
    hold on;plot(block{2}.t,block{2}.B1amp*1.5,'y','linewidth',2);
    hold off;
    axis([0,t(end),-0.2 1.5]);
    set(axes1, 'color', 'k');  
    set(axes1,'XTick',linspace(0,t(end),3),'FontSize',12,'FontWeight','bold','XColor','w');
    set(axes1,'xticklabels',{num2str(0),[num2str(t(end)/2*1e3),' ms'], [num2str(t(end)*1e3),' ms']});
    set(axes1,'YTick',[-0.2 0  0.2 0.4 0.6 0.8 1],'FontSize',12,'FontWeight','bold','YColor','w');
    text(t(end)/4,1.4, ['SINC pulse TA (', num2str(m),')'], 'FontSize',10,'color','r')
    text(t(end)/4,1.25, ['SLR pulse TA (', num2str(m),')'], 'FontSize',10,'color','y')

    axes2=subplot(122);
    plot(f,abs(spins1.M(:,1)+1i*spins1.M(:,2)),'r','linewidth',2);
    hold on;plot(f,abs(spins2.M(:,1)+1i*spins2.M(:,2)),'y','linewidth',2);

%     plot(f,abs(spins1.M(:,1)),'r','linewidth',2);
%     hold on;plot(f,abs(spins2.M(:,1)+1i*spins2.M(:,2)),'y','linewidth',2);
    hold off;

    % plot(f,spins1.M(:,3),'r','linewidth',2);
    % hold on;plot(f,spins2.M(:,3),'y','linewidth',2);

    f_bandwith=f(end)-f(1);

    text(-f_bandwith/4,1.4, 'SINC pulse excited profile', 'FontSize',10,'color','r')
    text(-f_bandwith/4,1.25, 'SLR pulse excited profile', 'FontSize',10,'color','y')

    axis([-f_bandwith/2 f_bandwith/2 -0.2,1.5]);
    %     title('Excited Magnetization(Mxy)','color','y');
    set(axes2, 'color', 'k');  
    set(axes2,'XTick',linspace(-f_bandwith/2,f_bandwith/2,3),'FontSize',12,'FontWeight','bold','XColor','w');
    set(axes2,'xticklabels',{num2str(-round(f_bandwith/2*1e-2)/10),'0', [num2str(round(f_bandwith/2*1e-2)/10),' kHz']});
    set(axes2,'FontSize',12,'FontWeight','bold','YColor','w');

    frame=getframe(gcf);
    writeVideo(writerObj,frame);
end        

close(writerObj);


