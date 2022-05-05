function pulse = PulseCreateRosenfeld2(TBNarrow, IB, omega0, numSteps)
% Creates an asymmetric excitation pulse based on 
%
%   Rosenfeld, Panfil and Zur, MRM 37:793-801 (1997): Design of adiabatic
%   pulses for fat-suppression using analytic solutions of the Bloch 
%   equation
%
% Input Variables
% Name       Description
% TBNarrow   Transition bandwidth, narrow side (Hz)
% TBWide     Transition bandwidth, wide side (Hz)
% IB         Inversion bandwidth (Hz)
% numSteps   Total number of pulse steps
%
% Output Variables
% Name       Description
% pulse      The output pulse structure

rhoMin = 11.09; 

T = 1.3/(TBNarrow/2);
b = T/171.33;
a = rhoMin*b;

f0 = 0;
S = [IB*T/2.255 + 7.625;
     f0*T/1.13 - 6.46];
M = [-1, 11.09;
     -1, -13.09];
x = M^(-1)*S;    
c = x(1);
d = x(2);

pulse = PulseCreateRosenfeld(a,b,c,d, omega0, numSteps);
