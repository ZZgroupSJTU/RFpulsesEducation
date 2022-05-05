function block=combineblocks(vargarin)


t=[];
B1amp=[];
B1phase=[];
Gx=[];
Gy=[];
Gz=[];
ADCflag=[];
for k=1:length(vargarin)
    if isempty(t)
        t=cat(2,t,vargarin{k}.t);
    else
        t=cat(2,t,t(end)+vargarin{k}.t);
    end
    B1amp=cat(2,B1amp,vargarin{k}.B1amp);
    B1phase=cat(2,B1phase,vargarin{k}.B1phase);
    Gx=cat(2,Gx,vargarin{k}.Gx);
    Gy=cat(2,Gy,vargarin{k}.Gy);
    Gz=cat(2,Gz,vargarin{k}.Gz);
%     ADCflag=cat(2,ADCflag,vargarin{k}.ADCflag);
end

block.tof=vargarin{1}.tof;
block.t=t;
block.B1amp=B1amp;
block.B1phase=B1phase;
block.Gx=Gx;
block.Gy=Gy;
block.Gz=Gz;
% block.ADCflag=all(ADCflag);
