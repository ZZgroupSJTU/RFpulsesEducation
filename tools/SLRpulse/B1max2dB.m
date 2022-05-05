function [PowerCoarse,PowerFine]=B1max2dB(B1max,varargin)

if isempty(varargin) || length(varargin)==1
    error('The reference powercal at least 2 input');
elseif length(varargin)==2
    refB1max=varargin{1};
    refpwr=varargin{2};
elseif length(varargin)>=3
    refB1max=varargin{1};
    refpwr=varargin{2}*varargin{3}/4095;
end

if refpwr>60
    warning(['The reference power ',num2str(refpwr),' is too large']);
end

if refB1max<1e-2
    warning('The refB1max is input as pw way');
    refB1max=1e-3/(4*refB1max); 
end

disp(['refB1max = ',num2str(refB1max), '   refpwr = ',num2str(refpwr)]);

rfpwr= refpwr+20*log10(B1max/refB1max);
rfpwr1=4095*10^(-(ceil(rfpwr)-rfpwr)/20);
PowerCoarse=ceil(rfpwr);
PowerFine=round(rfpwr1);

disp(['rfpwr = ',num2str(rfpwr)]);

if rfpwr>60
    warning(['The caculate power ',num2str(rfpwr),' is too large']);
end