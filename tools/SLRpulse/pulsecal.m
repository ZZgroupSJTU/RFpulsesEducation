function [Corsepower,Finepower,Power]=pulsecal(Refname,Refduration,Refpower,Refflip,Pulsename,Pulseduration,Pulseflip,varargin)


path='/home/Zhiyongi/vnmrsys/shapelib';
if ~isempty(varargin)
    path=varargin{1};
end


fulfilename=[path,'/',Refname,'.RF'];
fid = fopen(fulfilename,'r'); 
if fid==-1
    error('Error opening the aglient file');
end
Refintegral=freadpar(fid,'INTEGRAL',1);
fclose(fid);

RefB1max=10^3*(Refflip*pi/180)/Refintegral/Refduration/2/pi;       % kHz

fulfilename=[path,'/',Pulsename,'.RF'];
fid = fopen(fulfilename,'r'); 
if fid==-1
    error('Error opening the aglient file');
end
Pulseintegral=freadpar(fid,'INTEGRAL',1);
fclose(fid);

PulseB1max=10^3*(Pulseflip*pi/180)/Pulseintegral/Pulseduration/2/pi;       % kHz

Power=Refpower+20*log10(PulseB1max/RefB1max);

Corsepower=ceil(Power);
Finepower=4095*10^(-(Corsepower-Power)/20);
