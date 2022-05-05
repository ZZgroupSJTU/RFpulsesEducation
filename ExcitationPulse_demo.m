%% 
% This is a RF pulse simulation demo using a Bloch simulator
% Created for ISMRM enducational talk of RF pulse in May 2022
% Zhiyong Zhang @ zhiyong.zhang@sjtu.edu.cn
%%
clc;clear;close all;
addpath(genpath(pwd));
%% definite a RF profile information and create a pulse structure using Assaf's SpinTool
%  PulseCreateSinc create a sinc pulse
pulse_duration=0.0025;   %[ms]
pulse_number_of_point=128;
pulse_number_of_lobe=3;
pulse_bandwidth=2*pulse_number_of_lobe/pulse_duration;   %[Hz]
pulse_flipangle=90;
pulse_initphase=0;
pulse=PulseCreateSinc(pulse_number_of_lobe,...
                        pulse_duration*1e3,...
                        pulse_number_of_point,...
                        pulse_flipangle,...
                        pulse_initphase);
%% Simulate the spin dynamics without Gradient
% Defined a 1D spin system
f=[-2000,0,1000];                               % 3 resonace frequency [Hz]
spins0.M0=ones(length(f),1);                                            
spins0=Init_Spinsystem({'cs',f},spins0);

% create a sequence block with RF without a slice-selection gradient
block{1}=Generate_seqblock(0,...                                   %Reference frequency
                           pulse_duration,...                      %block's duration
                           pulse_number_of_point,...               %block's number of point for simulation
                           pulse.RFamp.*exp(1i*pulse.RFphase),...  % B1 field (complex value)
                           0,...                                   % Gx gradient (support constant value and a shape)
                           0,...                                   % Gy gradient (support constant value and a shape)
                           0);                                     % Gz gradient (support constant value and a shape)
% Show the dynamics
[spins1,Mevolution1]=ConsoleSimulator(spins0,block{1});
params.BGColor='k';
params.viewSeqblock=true;
params.SaveFormat='mp4';
params.filename='Excitation_noGz_demo';
figure();showMagnetizationWithSeqblock_Dynamic(Mevolution1.Mt(:,:,:),block{1},params)
%% Simulate the spin dynamics with slice-selection Gradient
% Let's re-initialize the spin system with z positions (slice dimension)
clear spins0;
Nz=256;
Lz=2;                                     %[cm];
z=linspace(-Lz/2,Lz/2-Lz/Nz,Nz);
spins0.M0=ones(Nz,1);                              
spins0.T1=1000*ones(Nz,1);                % Let's ignore the T1 relaxation for procession demo
spins0.T2=1000*ones(Nz,1);                % Let's ignore the T2 relaxation for procession demo
spins0=Init_Spinsystem({'z',z},spins0);
% create a sequence block with RF with a slice-selection gradient
thickness= 0.5;                            % [cm]
gammaHz=4257.4;                            % [Hz/guass]
Gz=pulse_bandwidth/gammaHz/thickness;      % calculate the slice selection gradient (gauss/cm)
block{1}=Generate_seqblock(0,...                                   %Reference frequency
                           pulse_duration,...                      %block's duration
                           pulse_number_of_point,...               %block's number of point for simulation
                           pulse.RFamp.*exp(1i*pulse.RFphase),...  % B1 field (complex value)
                           0,...                                   % Gx gradient (support constant value and a shape)
                           0,...                                   % Gy gradient (support constant value and a shape)
                           Gz);                                     % Gz gradient (support constant value and a shape)

% Show the excited slice profile
[spins2,Mevolution2]=ConsoleSimulator(spins0,block{1});
figure()
subplot(121);showseqblock(block(1),1);
subplot(122);plot(z(:),abs(spins2.M(:,1)+1i*spins2.M(:,2)));
% show the dynamics
params.BGColor='k';
params.viewSeqblock=true;
params.SaveFormat='mp4';
params.filename='Excitation_withGz_demo';
figure();showMagnetizationWithSeqblock_Dynamic(Mevolution2.Mt(1:4:end,:,1:1:end),block{1},params)
