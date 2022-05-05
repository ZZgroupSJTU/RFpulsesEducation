function pulse=pinsPulse(pulse,Nsubs,dc,maxslewrate,glim)
%%
% maxslewrate: [mT/m/ms]
% glim:  [mT/m]
pts=length(pulse.RFamp);
if (mod(pts,Nsubs)~=0)
    error('please create multiple-Nsubs pulse');
end
pindex=reshape(zeros(1,pts),[],Nsubs);
M=size(pindex,1);
pindex(round(M*dc):end,:)=1;
indexlogical=logical(reshape(pindex,1,[]));

RFInt0=sum(pulse.RFamp); 
pulse.RFamp(indexlogical)=0;
pulse.RFphase(indexlogical)=0;
RFInt1=sum(pulse.RFamp); 
pulse.RFamp=pulse.RFamp*RFInt0/RFInt1;
if nargin==3
    GxInt0=sum(pulse.Gx); 
    GyInt0=sum(pulse.Gy);
    GzInt0=sum(pulse.Gz);
    pulse.Gx(~indexlogical)=0;
    pulse.Gy(~indexlogical)=0;
    pulse.Gz(~indexlogical)=0;
    GxInt1=sum(pulse.Gx); 
    GyInt1=sum(pulse.Gy);
    GzInt1=sum(pulse.Gz);
    if GxInt1==0
        GxInt1=1;
    end
    if GyInt1==0
        GyInt1=1;
    end
    if GzInt1==0
        GzInt1=1;
    end
    pulse.Gx=pulse.Gx*GxInt0/GxInt1;
    pulse.Gy=pulse.Gy*GyInt0/GyInt1;
    pulse.Gz=pulse.Gz*GzInt0/GzInt1;

else
    dt= pulse.tp/length(pulse.RFamp);
    G=[pulse.Gx;pulse.Gy;pulse.Gz];
    Tsub=pulse.tp/Nsubs;
    bliparea=max(G,[],2)*Tsub;
    if length(maxslewrate)==1
        maxslewrate=maxslewrate*ones(1,3);
    elseif length(maxslewrate)>3
        maxslewrate=maxslewrate(1:3);
    end
    if length(glim)==1
        glim=glim*ones(1,3);
    elseif length(glim)>3
        glim=glim(1:3);
    end
    
    glim=glim/10;  % mT/m  --> G/cm
    maxslewrate=maxslewrate/10;      % mT/m/ms  ----> G/cm/ms;
    
    blipblock=zeros(1,M);
    for k=1:3
        if bliparea(k)>0
            %gflat*Tsub+0.5*Tsub*(1-dc)*(0.5*Tsub*(1-dc)*maxslewrate)=bliparea;
            gflat=(bliparea(k)-0.5*Tsub*(1-dc)*(0.5*Tsub*(1-dc)*maxslewrate(k)))/Tsub; 
            gmax=gflat+Tsub*(1-dc)/2*maxslewrate(k);
            if gflat<0
               gflat=0;
               %0.5*Tsub*(1-dc)*gmax=bliparea
               gmax=bliparea(k)/(0.5*Tsub*(1-dc));
            end
            if gmax>glim(k)
                gmax=glim(k);
                %gflat*Tsub+0.5*Tsub*(1-dc)*(gmax-gflat)=bliparea(k);
                gflat=(bliparea(k)-0.5*Tsub*(1-dc)*gmax)/(Tsub-0.5*Tsub*(1-dc));
            end
            temp=linspace(gflat,gmax,round(Tsub*(1-dc)/dt/2));
            temp=[temp,temp((end-1):-1:1)];
            blipblock(end:-1:(end-length(temp)+1))=temp;
            blipblock(1:(end-length(temp)))=gflat;
            G(k,:)=reshape(repmat(blipblock.',[1,Nsubs]),1,[]);
        end
    end
    pulse.Gx=G(1,:);
    pulse.Gy=G(2,:);
    pulse.Gz=G(3,:);
end