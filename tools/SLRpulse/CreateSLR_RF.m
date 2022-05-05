function p1_rf=CreateSLR_RF(SLR_RF,varargin)
%%
% SLRpulsetype='alpha' ,  'excitation' ,  'inversion',   'refocusing'
% SLRfiltertype= 'leastsquare' , 'equiripple' ,  'complexequiripple',  'quadraticequiripple' 
%%


pulname=SLR_RF.pulname;
duration=SLR_RF.duration;
pts=SLR_RF.pts;
flip=SLR_RF.flip;
bandwidth=SLR_RF.bandwith;
inripple=SLR_RF.inripple;
outripple=SLR_RF.outripple;
SLRfiltertype=SLR_RF.filtertype;
SLRpulsetype=SLR_RF.pulsetype;


refpw=10e-6;
reftpwr=58;
pulsepath=pwd;
if length(varargin)==1
    refpw=varargin{1};
elseif length(varargin)==2
    refpw=varargin{1};
    reftpwr=varargin{2};
elseif length(varargin)>=2
    refpw=varargin{1};
    reftpwr=varargin{2};
    pulsepath=varargin{3};
end



%%

p1_rf=init_rf(pulname);
p1_rf.duration=duration;
p1_rf.pts=pts;
p1_rf.flip=flip;
p1_rf.bandwidth=bandwidth;                          % in [kHz] : RF pulse bandwidth

%%
if isfield(SLR_RF,'multiband') && strcmpi(SLR_RF.multiband,'true')
    % multif:  frequency
    % multiBW: BandWith
    % multiTW: Transition bands
    
    f=2*reshape(SLR_RF.multif,1,[]);
    BW=reshape(SLR_RF.multiBW,1,[]);
    TW=reshape(SLR_RF.multiTW,1,[]);
    SLRParams.filtern = pts-1;
    SLRParams.polyn = pts;
    SLRParams.bandf = [-1 reshape([f-BW/2-TW/2;f-BW/2+TW/2;f+BW/2-TW/2;f+BW/2+TW/2]*duration/pts,1,[]) 1];
    SLRParams.banda = [0 reshape([0.0001*ones(1,length(f));ones(1,length(f))*sin(pi*flip/360);...
                       ones(1,length(f))*sin(pi*flip/360);0.0001*ones(1,length(f))],1,[]),0];
    SLRParams.bandw = [SLR_RF.multiws,reshape([ones(1,length(f))*SLR_RF.multiwp; ones(1,length(f))*SLR_RF.multiws],1,[])];  
    SLRParams.filtertype='complexequiripple';
else
    SLRParams=init_SLRparam(p1_rf.duration,p1_rf.bandwidth,inripple,outripple,SLRfiltertype,SLRpulsetype,p1_rf.pts,p1_rf.flip);
end
    
[a, b] = calc_SLRpolyn(SLRParams);

[thetalong, philong] = inverseSLR(a,b);

p1_rf=calc_SLRrf(p1_rf,thetalong,philong,refpw,reftpwr);

% saveRFfile(p1_rf,pulsepath);