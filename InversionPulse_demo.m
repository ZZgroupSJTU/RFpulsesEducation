%% 
% This is a RF pulse simulation demo using a Bloch simulator
% Created for ISMRM enducational talk of RF pulse in May 2022
% Zhiyong Zhang @ zhiyong.zhang@sjtu.edu.cn
%% Defined a 1D spin system
clc;clear;close all;
%% initialize the spin system
f=[0,0];                                               % Two water spins with same frequency：one for regular tissue, one for fluid
spins0.T1=[0.8;2.8];                                     % T1=500 ms for regular tissue, T1=2000 ms for fluid
spins0.M0=ones(length(f),1);
spins0=Init_Spinsystem({'cs',f},spins0);
%% definite a RF profile information and create a pulse structure using Assaf's SpinTool
pulse_duration=0.0025;   %[ms]
pulse_number_of_point=256;
pulse_number_of_lobe=3;
pulse_bandwidth=2*pulse_number_of_lobe/pulse_duration;   %[Hz]
pulse_flipangle=180;
pulse_initphase=0;
pulse=PulseCreateSinc(pulse_number_of_lobe,...
                        pulse_duration*1e3,...
                        pulse_number_of_point,...
                        pulse_flipangle,...
                        pulse_initphase);
%% create a sequence block with RF without a slice-selection gradient
block{1}=Generate_seqblock(0,...                                 %Reference frequency
                           pulse_duration,...                    %block's duration
                           pulse_number_of_point,...             %block's number of point for simulation
                           pulse.RFamp.*exp(1i*pulse.RFphase),...  % B1 field (complex value)
                           0,...                                 % Gx gradient (support constant value and a shape)
                           0,...                                 % Gy gradient (support constant value and a shape)
                           0);                                   % Gz gradient (support constant value and a shape)

TI=2;
block{2}=Generate_seqblock(0,...                                 %Reference frequency
                           TI,...                                %block's duration
                           128,...                                %block's number of point for simulation
                           0,...                                 % B1 field (complex value)
                           0,...                                 % Gx gradient (support constant value and a shape)
                           0,...                                 % Gy gradient (support constant value and a shape)
                           0);                                   % Gz gradient (support constant value and a shape)

                      
block{3}=Generate_seqblock(0,...                                 %Reference frequency
                           pulse_duration,...                    %block's duration
                           pulse_number_of_point,...             %block's number of point for simulation
                           pulse.RFamp/2.*exp(1i*pulse.RFphase),...  % B1 field (complex value)
                           0,...                                 % Gx gradient (support constant value and a shape)
                           0,...                                 % Gy gradient (support constant value and a shape)
                           0);                                   % Gz gradient (support constant value and a shape)
seqblock=combineblocks(block);
%% Show the first inversion pulse dynamics
spins0.M=repmat([0 0 1],[length(f),1]);
[spins1,Mevolution1]=ConsoleSimulator(spins0,block{1});
params.BGColor='k';
params.viewSeqblock=true;
params.SaveFormat='mp4';
params.filename='Inversion_pulse_demo';
figure();showMagnetizationWithSeqblock_Dynamic(Mevolution1.Mt(:,:,:),block{1},params)
%% Show the recovery dynamics
spins0.M=repmat([0 0 -1],[length(f),1]);
[spins2,Mevolution2]=ConsoleSimulator(spins0,block{2});
params.BGColor='k';
params.viewSeqblock=true;
params.SaveFormat='Inversion_recovery_demo';
figure();showMagnetizationWithSeqblock_Dynamic(Mevolution2.Mt(:,:,:),block{2},params)
%% Show the inversion puslse+TI+ exciation puse
spins0.M=repmat([0 0 1],[length(f),1]);
[spins3,Mevolution3]=ConsoleSimulator(spins0,seqblock);
params.BGColor='k';
params.viewSeqblock=true;
params.SaveFormat='mp4';
params.filename='Inversion_TI_M_demo';
figure();showMagnetizationWithSeqblock_Dynamic(Mevolution3.Mt(:,:,:),seqblock,params)
%%