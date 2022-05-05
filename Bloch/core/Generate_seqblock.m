function seqblock=Generate_seqblock(tof,Tb,Nb,rf,Gxgrad,Gygrad,Gzgrad)
%%  seqblock defined as: 
%                         seqblock.tof          [Hz]
%                         seqblock.t:           [s]
%                         seqblock.B1amp:       [kHz]
%                         seqblock.B1phase:     [rad]
%                         seqblock.Gx           [Gauss/cm]
%                         seqblock.Gy           [Gauss/cm]
%                         seqblock.Gz           [Gauss/cm]
% Zhiyong Zhang@June 1,2016, zhiyongxmu@gmail.com; zhiyong.zhang@weizmann.ac.il   
%%
seqblock.tof=tof;
seqblock.t=linspace(0,Tb,Nb);
if nargin<3
    error('~come on, at least give me three inputs');
end

seqblock.B1amp=zeros(1,Nb);
seqblock.B1phase=zeros(1,Nb);
seqblock.Gx=zeros(1,Nb);
seqblock.Gy=zeros(1,Nb);
seqblock.Gz=zeros(1,Nb);

if nargin>3;
   if ~isstruct(rf)
       if numel(rf)==1
           seqblock.B1amp=ones(1,Nb)*abs(rf);
           seqblock.B1phase=ones(1,Nb)*angle(rf);
       else
           seqblock.B1amp=interp1(linspace(0,Tb,length(rf)),abs(rf),seqblock.t);
           seqblock.B1phase=interp1(linspace(0,Tb,length(rf)),angle(rf),seqblock.t);
       end
   else
        if rf.duration > Tb
            error('Block duration is less the duration of RF pulse');
        end
        if rf.duration == 0
           rf.duration=Tb;
        end 
        B1=rf.B1max.*rf.B1.*exp(1i*2*pi*rf.channelphase/360);
        if (abs(rf.B1(1))==1)
            B1(1)=0;
        end
        nfix=round(length(B1)*(Tb-rf.duration)/rf.duration/2);
        B1=[zeros(1,nfix),B1,zeros(1,nfix)];
        B1=interp1(linspace(0,Tb,length(B1)),B1,seqblock.t);
        seqblock.B1amp=abs(B1);
        seqblock.B1phase=angle(B1);
   end
end
if nargin>4;
   if ~isstruct(Gxgrad)
       if numel(Gxgrad)==1
           seqblock.Gx=ones(1,Nb)*Gxgrad;
       else
           seqblock.Gx=interp1(linspace(0,Tb,length(Gxgrad)),Gxgrad,seqblock.t);
       end
   else
        if Gxgrad.duration > Tb
            error('Block duration is less the duration of Gxgrad');
        end
        if Gxgrad.duration == 0
           Gxgrad.duration=Tb;
        end 
        Gx=Gxgrad.ampdata.*Gxgrad.head.strength/max(Gxgrad.ampdata);
        nfix=round(length(Gx)*(Tb-Gxgrad.duration)/Gxgrad.duration/2);
        Gx=[zeros(1,nfix),Gx,zeros(1,nfix)];
        seqblock.Gx=interp1(linspace(0,Tb,length(Gx)),Gx,seqblock.t);
        if isfield(Gxgrad,'negative') && Gxgrad.negative
           seqblock.Gx=-seqblock.Gx;
        end
   end
end

if nargin>5;
   if ~isstruct(Gygrad)
       if numel(Gygrad)==1
           seqblock.Gy=ones(1,Nb)*Gygrad;
       else
           seqblock.Gy=interp1(linspace(0,Tb,length(Gygrad)),Gygrad,seqblock.t);
       end
   else
        if Gygrad.duration > Tb
            error('Block duration is less the duration of Gygrad');
        end
        if Gygrad.duration == 0
           Gygrad.duration=Tb;
        end 
        Gy=Gygrad.ampdata.*Gygrad.head.strength/max(Gygrad.ampdata);
        nfix=round(length(Gy)*(Tb-Gygrad.duration)/Gygrad.duration/2);
        Gy=[zeros(1,nfix),Gy,zeros(1,nfix)];
        seqblock.Gy=interp1(linspace(0,Tb,length(Gy)),Gy,seqblock.t);
        if isfield(Gygrad,'negative') && Gygrad.negative
           seqblock.Gy=-seqblock.Gy;
        end
   end
end

if nargin>6;
   if ~isstruct(Gzgrad)
       if numel(Gzgrad)==1
           seqblock.Gz=ones(1,Nb)*Gzgrad;
       else
           seqblock.Gz=interp1(linspace(0,Tb,length(Gzgrad)),Gzgrad,seqblock.t);
       end
   else
        if Gzgrad.duration > Tb
            error('Block duration is less the duration of Gygrad');
        end
        if Gzgrad.duration == 0
           Gzgrad.duration=Tb;
        end 
        Gz=Gzgrad.ampdata.*Gzgrad.head.strength/max(Gzgrad.ampdata);
        nfix=round(length(Gz)*(Tb-Gzgrad.duration)/Gzgrad.duration/2);
        Gz=[zeros(1,nfix),Gz,zeros(1,nfix)];
        seqblock.Gz=interp1(linspace(0,Tb,length(Gz)),Gz,seqblock.t);
        if isfield(Gzgrad,'negative') &&Gzgrad.negative
           seqblock.Gz=-seqblock.Gz;
        end
   end
end

