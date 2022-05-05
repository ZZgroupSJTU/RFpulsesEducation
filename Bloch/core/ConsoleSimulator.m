function [spinsout,Mevolution]=ConsoleSimulator(spinsin,seqblock)

%%  input: 
%             seqblock.tof          [Hz]
%             seqblock.t:           [s]
%             seqblock.B1amp:       [kHz]
%             seqblock.B1phase:     [rad]
%             seqblock.Gx           [Gauss/cm]
%             seqblock.Gy           [Gauss/cm]
%             seqblock.Gz           [Gauss/cm]
%  Output:
%             spins.rcs   (nspins*4)    [cm,cm,cm,Hz]
%             spins.M0    (nspins*3)    
%             spins.M     (nspins*3)
%             spins.T1    (nspins*1)    [s]
%             spins.T2    (nspins*1)    [s]
%             spins.B0    (nspins*1)    [Hz]
%             spins.B1    (nspins*1)    [Hz]
%             spins.gammaHz             [Hz/gauss]
%
%             Mevolution.t=t;
%             Mevolution.Mxyt   nspins*Nt;  
%             Mevolution.Mzt    nspins*Nt;   
% Zhiyong Zhang@June 1,2016, zhiyongxmu@gmail.com; zhiyong.zhang@sjtu.edu.cn   
%%
if ~isfield(spinsin,'gammaHz')
    spinsin.gammaHz=4257.4;
end
gammaHz=spinsin.gammaHz;

t=seqblock.t;
Nt=length(t);
B1amp=seqblock.B1amp;
B1phase=seqblock.B1phase;

G=[gammaHz*seqblock.Gx;gammaHz*seqblock.Gy;gammaHz*seqblock.Gz;ones(1,Nt)];
nspins=size(spinsin.M0,1);

M0=spinsin.M0;
M=spinsin.M;
T1=spinsin.T1;
T2=spinsin.T2;

% B0=repmat(spinsin.B0,[1,length(t)])+spinsin.rcs*G+seqblock.tof;  %[Hz];
% evoangle=-2*pi*B0*tstep;

% cosevoangle=cos(evoangle);
% sinevoangle=sin(evoangle);


if nargout==2 
   Mevolution.t=t;
   Mevolution.B0t=zeros(nspins,Nt);
   Mevolution.Mt=zeros(nspins,3,Nt);
   Mevolution.Mxyt=complex(zeros(1,Nt));
   Mevolution.Mzt=complex(zeros(1,Nt));
   Mevolution.Mt(:,:,1)=M;
   Mevolution.Mxyt(:,1)=sum(M(:,1)+1i*M(:,2));
   Mevolution.Mzt(:,:,1)=sum(M(:,3));
end
for m=2:Nt 
%     tic;
    tstep=t(m)-t(m-1);
    smallflip=2*pi*B1amp*tstep*1e3;
    
    T2decay=exp(-tstep./T2);
    T1decay=exp(-tstep./T1);

    cosB1phase=cos(B1phase);
    sinB1phase=sin(B1phase);
    cossmallflip=cos(smallflip);
    sinsmallflip=sin(smallflip);

    R=zeros(3,3,Nt);
    R(1,1,:)=cosB1phase.^2+cossmallflip.*sinB1phase.^2;
    R(2,2,:)=cossmallflip.*cosB1phase.^2+sinB1phase.^2;
    R(3,3,:)=cossmallflip;
    R(1,2,:)=cosB1phase.*sinB1phase.*(cossmallflip-1);
    R(2,1,:)=R(1,2,:);
    R(1,3,:)=-sinB1phase.*sinsmallflip;
    R(3,1,:)=-R(1,3,:);
    R(2,3,:)=-cosB1phase.*sinsmallflip;
    R(3,2,:)=-R(2,3,:);

    M=M*R(:,:,m);

    B0tmp=spinsin.B0+spinsin.rcs*G(:,m)+seqblock.tof;
    Mx=(M(:,1).*cos(-2*pi*B0tmp*tstep)-M(:,2).*sin(-2*pi*B0tmp*tstep)).*T2decay;
    My=(M(:,1).*sin(-2*pi*B0tmp*tstep)+M(:,2).*cos(-2*pi*B0tmp*tstep)).*T2decay;

    Mz=M0(:,3)+(M(:,3)-M0(:,3)).*T1decay;
    M=[Mx,My,Mz];
    if nargout==2 
       Mevolution.B0t(:,m)=B0tmp; 
       Mevolution.Mt(:,:,m)=M; 
       Mevolution.Mxyt(:,m)=sum(Mx(:)+1i*My(:)); 
       Mevolution.Mzt(:,m)=sum(Mz(:)); 
    end
%     disp([num2str(m) ':' num2str(toc)]);
end
spinsout=spinsin;
spinsout.M=M;

