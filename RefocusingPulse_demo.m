%% 
% This is a RF pulse simulation demo using a Bloch simulator
% Created for ISMRM enducational talk of RF pulse in May 2022
% Zhiyong Zhang @ zhiyong.zhang@sjtu.edu.cn
%% Defined a 1D spin system
clc;clear;close all;
%% Use 11 spins to demostrate the backgroud field evolution of spins
nspins=11;
f=zeros(nspins,1);            
spins0.M0=ones(nspins,1);
spins0.B0=linspace(-10,10,nspins).';
spins0.T1=10000;
spins0.T2=10000;
spins0=Init_Spinsystem({'cs',f},spins0);
%% definite a RF profile information and create a pulse structure using Assaf's SpinTool
pulse_duration=0.004;   %[ms]
pulse_number_of_point=100;
pulse_number_of_lobe=3;
pulse_flipangle=90;
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

TEhalf=0.020;
block{2}=Generate_seqblock(0,...                                 %Reference frequency
                           TEhalf,...                                %block's duration
                           64,...                                %block's number of point for simulation
                           0,...                                 % B1 field (complex value)
                           0,...                                 % Gx gradient (support constant value and a shape)
                           0,...                                 % Gy gradient (support constant value and a shape)
                           0);                                   % Gz gradient (support constant value and a shape)

                      
block{3}=Generate_seqblock(0,...                                 %Reference frequency
                           pulse_duration,...                    %block's duration
                           pulse_number_of_point,...             %block's number of point for simulation
                           2*pulse.RFamp.*exp(1i*pulse.RFphase),...  % B1 field (complex value)
                           0,...                                 % Gx gradient (support constant value and a shape)
                           0,...                                 % Gy gradient (support constant value and a shape)
                           0);                                   % Gz gradient (support constant value and a shape)

block{4}=Generate_seqblock(0,...                                 %Reference frequency
                           TEhalf+pulse_duration/2,...                                %block's duration
                           64,...                                %block's number of point for simulation
                           0,...                                 % B1 field (complex value)
                           0,...                                 % Gx gradient (support constant value and a shape)
                           0,...                                 % Gy gradient (support constant value and a shape)
                           0);                                   % Gz gradient (support constant value and a shape)     
seqblock=combineblocks(block);
%%
spins0.M=repmat([0 0 1],[length(f),1]);
[spins,Mevolution]=ConsoleSimulator(spins0,seqblock);

params.BGColor='k';
params.viewSeqblock=true;
params.viewTrace=true;
% params.SaveFormat='mp4';
params.SaveFormat='none';
params.filename='Refocusing_pulse_90x180x_demo';
figure();showMagnetizationWithSeqblock_Dynamic(Mevolution.Mt(:,:,:),seqblock,params)
%% 90x+ 180y

block{3}=Generate_seqblock(0,...                                 %Reference frequency
                           pulse_duration,...                    %block's duration
                           pulse_number_of_point,...             %block's number of point for simulation
                           2*pulse.RFamp.*exp(1i*(pulse.RFphase+pi/2)),...  % B1 field (complex value)
                           0,...                                 % Gx gradient (support constant value and a shape)
                           0,...                                 % Gy gradient (support constant value and a shape)
                           0);                                   % Gz gradient (support constant value and a shape)
  
seqblock=combineblocks(block);

spins0.M=repmat([0 0 1],[length(f),1]);
[spins,Mevolution]=ConsoleSimulator(spins0,seqblock);

params.BGColor='k';
params.viewSeqblock=true;
params.viewTrace=true;
% params.SaveFormat='mp4';
params.SaveFormat='none';
params.filename='Refocusing_pulse_90x180y_demo';
figure();showMagnetizationWithSeqblock_Dynamic(Mevolution.Mt(:,:,:),seqblock,params)