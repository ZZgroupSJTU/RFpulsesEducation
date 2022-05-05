function spins = Init_Spinsystem(Cellinput,spins)
%%
%  Initilize the spin system for simulation
%  The input dim is a cell as following:
%        {'z',zpostion}
%        {'z',zposition,'x',xposition};
%        {'z',zposition,'x',xposition,'y',yposition};
%        {'z',zposition,'x',xposition,'y',yposition,'cs',foffset};
%  Output spins: 
%             spins.rcs   (nspins*4)    [cm,cm,cm,Hz]
%             spins.M0    (nspins*3)    
%             spins.M     (nspins*3)
%             spins.T1    (nspins*1)    [s]
%             spins.T2    (nspins*1)    [s]
%             spins.B0    (nspins*1)    [Hz]
%             spins.B1    (nspins*1)    [Hz]
%             spins.gammaHz             [Hz/gauss]
% Zhiyong Zhang@June 1,2016, zhiyongxmu@gmail.com; zhiyong.zhang@weizmann.ac.il        
%%

Ncell=length(Cellinput);
Nx=1;Ny=1;Nz=1;Ncs=1;
spins.gammaHz=4257.4;
for k=1:2:Ncell
    switch Cellinput{k}
        case 'x'
            xtemp(:,1,1,1)=Cellinput{k+1};
            Nx=length(xtemp);
        case 'y'
            ytemp(1,:,1,1)=Cellinput{k+1};
            Ny=length(ytemp);
        case 'z'
            ztemp(1,1,:,1)=Cellinput{k+1};
            Nz=length(ztemp);
        case 'cs'
            cstemp(1,1,1,:)=Cellinput{k+1};
            Ncs=length(cstemp);
    end
end
if ~exist('xtemp','var')
    xtemp=0;
end
if ~exist('ytemp','var')
    ytemp=0;
end
if ~exist('ztemp','var')
    ztemp=0;
end
if ~exist('cstemp','var')
    cstemp=0;
end
nspins=Nx*Ny*Nz*Ncs;
spins.rcs=zeros(nspins,4);
spins.rcs(:,1)=reshape(repmat(xtemp,1,Ny,Nz,Ncs),[],1);
spins.rcs(:,2)=reshape(repmat(ytemp,Nx,1,Nz,Ncs),[],1);
spins.rcs(:,3)=reshape(repmat(ztemp,Nx,Ny,1,Ncs),[],1);
spins.rcs(:,4)=reshape(repmat(cstemp,Nx,Ny,Nz,1),[],1);
           
if ~isfield(spins,'M0')
    spins.M0=repmat([0,0,1],nspins,1);
end
if ~isfield(spins,'T1')
    spins.T1=1*ones(nspins,1);
end
if ~isfield(spins,'T2')
    spins.T2=0.5*ones(nspins,1);
end
if ~isfield(spins,'B0')
    spins.B0=zeros(nspins,1);
end    
if ~isfield(spins,'B1')
    spins.B1=ones(nspins,1);
end


if size(spins.M0,2)==1
   spins.M0=[zeros(nspins,2),spins.M0];
end

if size(spins.T1,1)==1
   spins.T1=spins.T1*ones(nspins,1);
end

if size(spins.T2,1)==1
   spins.T2=spins.T2*ones(nspins,1);
end


if size(spins.M0,1)~=size(spins.rcs)
    error('please input same size MO with rcs');
end

if size(spins.T1,1)~=size(spins.rcs)
    error('please input same size T1 with rcs');
end

if size(spins.T2,1)~=size(spins.rcs)
    error('please input same size T2 with rcs');
end

if size(spins.B0,1)~=size(spins.rcs)
    error('please input same size B0 with rcs');
end

if size(spins.B1,1)~=size(spins.rcs)
    error('please input same size B1 with rcs');
end

spins.M=spins.M0;


% spins.rcs = gpuArray(spins.rcs );
% spins.M0 = gpuArray(spins.M0 );
% spins.M = gpuArray(spins.M );
% spins.B0 = gpuArray(spins.B0 );
% spins.B1 = gpuArray(spins.B1 );
% spins.T1 = gpuArray(spins.T1 );
% spins.T2 = gpuArray(spins.T2 );