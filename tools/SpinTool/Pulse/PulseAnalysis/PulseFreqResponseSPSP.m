function [vv,zz,Mxy]=PulseFreqResponseSPSP(pulse,sw,Nv,FOV,Nz,M0);
% SYNTAX: Mxy=PulseFreqResponseSPSP(pulse,sw,Nv,FOV,Nz,M0);
%
% Much like PulseFreqResponse, only this creates a 2D contour plot
% with one axis for frequency and the other for position (z).

% Step I: Create a spin structure of the following form:
% [zi .... zf,  zi .... zf,   zi .... zf, ............ ]
%    1st cs       2nd cs         3rd cs

vv = linspace(-sw/2,sw/2,Nv);
zz  = linspace(-FOV/2,FOV/2,Nz);
counter   = 0;
for k=1:Nv
    for p=1:Nz
        counter = counter + 1;
        spinstemp(counter).r = [0; 0; zz(p)];   % mm
        spinstemp(counter).M = M0;
        spinstemp(counter).cs = vv(k);   % kHz
    end;
end;

% Apply pulse
spinsout = ApplyPulse(spinstemp,pulse);

% Extract output magnetization
Mxy = zeros(Nv,Nz);
counter = 0;
for k=1:Nv
    for p=1:Nz
        counter = counter + 1;
        Mxy(k,p) = spinsout(counter).M(1) + i*spinsout(counter).M(2);
    end;
end;