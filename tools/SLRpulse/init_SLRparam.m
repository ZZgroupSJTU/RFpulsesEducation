function SLRParams=init_SLRparam(T,B,inripple,outripple,SLRfiltertype,SLRpulsetype,RF_nsteps,RFtheta)

% Calculate the filter design parameters from the pulse design parameters,
% using the Pauly relations updated by Le Roux in his thesis.


switch SLRpulsetype
    
    case 'alpha'
         
        passripple = inripple;
        stopripple = outripple;
        
    case 'excitation'
        
        passripple = 2*inripple;
        stopripple = outripple/sqrt(2);   
        
    case 'inversion'
        
        passripple = inripple/2;
        stopripple = sqrt(outripple);
        
    case 'refocusing'

        passripple = inripple/2;
        stopripple = sqrt(outripple);
        
    otherwise
        
        error('Unknown pulse type')
        
end

%--------------------calc f_p,f_s-------------------------%
%        B = (F_p+F_s)/delta_t         
%        ftw = (F_p-F_s)/(F_p+F_s) 
%---------------------------------------------------------%

Dinf=-0.6*log(sqrt(passripple*stopripple))-0.8;
ftw=Dinf/(T*B*1e3);                                              % fraction transition width
delta_t=T/(RF_nsteps-1);
f_p = 1e3*(B*delta_t-B*ftw*delta_t)/2;
f_s = 1e3*(B*delta_t+B*ftw*delta_t)/2;

w_p = 1;
% w_s = w_p*passripple/stopripple;
w_s = 1;
if(RFtheta>2*pi);
    RFtheta=RFtheta*pi/180;
end
sinhalftheta = sin(RFtheta/2);

% Prepare the input of the FIR filter design algorithms
SLRParams.filtern = RF_nsteps-1;
SLRParams.polyn = RF_nsteps;
SLRParams.bandf = [0 f_p f_s 1];
SLRParams.banda = [sinhalftheta sinhalftheta 0 0];
SLRParams.bandw = [w_p w_s];  
SLRParams.filtertype=SLRfiltertype;


