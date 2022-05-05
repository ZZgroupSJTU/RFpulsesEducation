function p2_rf=PolychromePulse(shape,duration,fcenter,flip,varargin)


pulsepath='/home/Zhiyongi/vnmrsys/shapelib';
if ~isempty(varargin)
    pulsepath=varargin{1};
end

[rfdata,p1_rf]=fread_agilentfile(pulsepath,shape,'RF');
p1_rf.duration=duration;
p1_rf=calc_rfpower(p1_rf,flip,duration);
p2_rf=p1_rf;
t=linspace(0,duration,length(rfdata(:,1))+1);
t=t(1:end-1);
rfdata=rfdata.';
B1sum=0;
for k=1:length(fcenter)
    B1sum=B1sum+rfdata(2,:).*exp(1i*2*pi*rfdata(1,:)/360).*exp(1i*2*pi*fcenter(k)*t);
end

p2_rf.B1max=p1_rf.B1max*max(abs(B1sum))/max(abs(rfdata(2,:)));


p2_rf.amp=abs(B1sum)*1023/max(abs(B1sum));
p2_rf.phase=mod(phase(B1sum)*180/pi,360);
p2_rf.pulseName=['Poly_',shape];
saveRFfile(p2_rf,pulsepath);