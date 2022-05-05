%% 
% This is a RF pulse simulation demo using a Bloch simulator
% Created for ISMRM enducational talk of RF pulse in May 2022
% Zhiyong Zhang @ zhiyong.zhang@sjtu.edu.cn
%%
clc;clear;close all;
addpath(genpath(pwd));
%% Defined a 1D spin system
Nz=256;                                   % 256 points along z axis
Lz=1;                                     % Unit: [cm];
z=linspace(-Lz/2,Lz/2-Lz/Nz,Nz);          % Physical positions
spins0.M0=ones(Nz,1);                     % This is generally to create the proton density contrast                                       
spins0.T1=1000*ones(Nz,1);                % Let's ignore the T1 relaxation for procession demo, default value is 1; unit [sec]
spins0.T2=1000*ones(Nz,1);                % Let's ignore the T2 relaxation for procession demo, default value is 0.5; unit [sec]
spins0=Init_Spinsystem({'z',z},spins0);   % initilize the spin system
%% Defined a evolution block without any RF and gradient, therefore, only precession is there
tevolution=10e-6;                                              %  [s]: evolution duration   
measure_point=400;                                             %  measure points
block{1}=Generate_seqblock(0,...                               %  Reference frequency 
                           tevolution,...                      %  block's duration
                           measure_point,...                   %  block's number of point for simulation
                           0,...                               %  B1 field (complex value)
                           0,...                               %  Gx gradient (support constant value and a shape)
                           0,...                               %  Gy gradient (support constant value and a shape)
                           0);                                 %  Gz gradient (support constant value and a shape)

%% Let look at the procession that the magnetizations are all along a vector (0  0.6 0.8);
spins0.M=repmat([0 0.6 0.8],[Nz,1]);                            % set the initial Magnetization to a vector (0 0.6 0.8) instead of in the equilibrium state (0 0 1)
[spins1,Mevolution1]=ConsoleSimulator(spins0,block{1});

params.BGColor='k';
params.viewSeqblock=false;
params.SaveFormat='mp4';
params.filename='precession_M';
figure();  
showMagnetizationWithSeqblock_Dynamic(Mevolution1.Mt(1,:,:),block{1},params)  %Noting only look at single spin magnetization
%% Defined a evolution block with a constant RF with 90 degree but without gradient
tevolution=6e-6;                                                %  [s]: evolution duration   
measure_point=200;                        
B1=0.25/tevolution*1e-3*exp(1i*2*pi*spins0.B0(1)*linspace(0,tevolution,measure_point));
block{1}=Generate_seqblock(0,...                               %  Reference frequency 
                           tevolution,...                      %  block's duration
                           measure_point,...                   %  block's number of point for simulation
                           B1,...                              %  B1 field (complex value)
                           0,...                               %  Gx gradient (support constant value and a shape)
                           0,...                               %  Gy gradient (support constant value and a shape)
                           0);                                 %  Gz gradient (support constant value and a shape)
%% Assume the larmor frequency of 1 MHz
spins0.B0=1e6*ones(Nz,1);                                      % Here, this is to simulate the larmor frequency (Assume w0=1 MHz) for this demo,
                                                               % but generally this is used to define the field inhomogeneity    
spins0.M=repmat([0 0 1],[Nz,1]);                               % reset the initial Magnetization to the equilibrium state 
[spins2,Mevolution2]=ConsoleSimulator(spins0,block{1});

params.BGColor='k';
params.viewSeqblock=false;
params.SaveFormat='mp4';
params.filename='Nutation_Mxyz_without_rotation frame';
figure();showMagnetizationWithSeqblock_Dynamic(Mevolution2.Mt(1,:,:),block{1},params); %Noting only look at single spin magnetization
%% In rotation frame with a frequency=larmor frequency
spins0.B0=0;                                                   % but generally this is used to define the field inhomogeneity    
spins0.M=repmat([0 0 1],[Nz,1]);                               % reset the initial Magnetization to the equilibrium state 

% B1 in rotation frame
tevolution=6e-6;                                                %  [s]: evolution duration   
measure_point=200;                        
B1=0.25/tevolution*1e-3*exp(1i*2*pi*spins0.B0(1)*linspace(0,tevolution,measure_point));
block{1}=Generate_seqblock(0,...                               %  Reference frequency 
                           tevolution,...                      %  block's duration
                           measure_point,...                   %  block's number of point for simulation
                           B1,...                              %  B1 field (complex value)
                           0,...                               %  Gx gradient (support constant value and a shape)
                           0,...                               %  Gy gradient (support constant value and a shape)
                           0);                                 %  Gz gradient (support constant value and a shape)
  
[spins3,Mevolution3]=ConsoleSimulator(spins0,block{1});

params.BGColor='k';
params.viewSeqblock=false;
params.SaveFormat='mp4';
params.filename='Nutation_Mxyz_in_rotation frame';
figure();showMagnetizationWithSeqblock_Dynamic(Mevolution3.Mt(1,:,:),block{1},params); %Noting only look at single spin magnetization

