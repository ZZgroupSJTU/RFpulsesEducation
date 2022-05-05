function pulse=pinsPulse(pulse,Nsubs,dc)
pts=length(pulse.RFamp);
if (mod(pts,Nsubs)~=0)
    error('please create multiple-Nsubs pulse');
end
pindex=reshape(zeros(1,pts),[],Nsubs);
M=size(pindex,1);
pindex(round(M*dc):end,:)=1;
indexlogical=logical(reshape(pindex,1,[]));

pulse.RFamp(indexlogical)=0;
pulse.RFphase(indexlogical)=0;
pulse.Gx(indexlogical)=0;
pulse.Gy(indexlogical)=0;
pulse.Gz(indexlogical)=0;