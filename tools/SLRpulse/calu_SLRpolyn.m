function [a,b] = calc_SLRpolyn(SLRParams,varargin)
% calu_SLRpolyn uses FIR filter design to find the B polynomial for a
% given target magnetisation and then obtains the corresponding
% minimum-phase A polynomial
% The SLR parameters are:
% filtertype: 'equiripple', 'complexequiripple' or 'leastsquare'
% filtern: order of the filter
% polyn: final size of the pulse, can be different from filtern for k not 0 
% bandf: limits of the band
% bandw: weight of the bands
% banda: amplitude of the bands
% k: amount of quadratic phase (not used for 'linear')

%%
% Load the relevant parameters
filtern = SLRParams.filtern;
bandf = SLRParams.bandf;
banda = SLRParams.banda;
bandw = SLRParams.bandw;
polyn = SLRParams.polyn;

% Create the B polynomial and calculate on the unit circle with high resolution

largen = max(20000,polyn);
switch SLRParams.filtertype

    case 'leastsquare'
        
        % Use the least-square algorithm 
        b = firls(filtern,bandf,banda,bandw);
    
    case 'equiripple'

        % Use the Parks-McClellan algorithm
        [b delta opt]= firpm(filtern,bandf,banda,bandw);
    
    case 'complexequiripple'

        % Use the complex Parks-McClellan algorithm
        [b delta opt] = ...
            cfirpm(filtern,bandf,{@multiband, banda},bandw);%,'skip_stage2');    
  
    case 'quadraticequiripple'    % Shulte method (2007)
        
        % Use the complex Parks-McClellan algorithm 
        % with a user-defined function.
        % The field k is removed to avoid applying the quadratic-phase
        % twice.
        if ~isempty(varargin)
            k = varargin{1};
            SLRParams.k=k;
            SLRParams = rmfield(SLRParams,'k');
        else 
            error('The quadratic phase k should be defined')
        end       
        [b delta opt] = ...
            cfirpm(filtern,bandf,{@quadraticband, banda, k},'skip_stage2');        
        
    otherwise
        
        error('unknown design method')
        
end

% Evalutate the B polynomial on the unit circle
% Three options are available to zero-fill the polynomial; the
% corresponding commands are found at the end of this function
blong = zeros(1,largen);
blong(largen/2-(filtern+1)/2+1:largen/2+(filtern+1)/2)=b; % 
Bunit = fftshift(fft(blong));


% Optionally add some quadratic phase
if isfield(SLRParams,'k')
    k = SLRParams.k;
    freqaxis = linspace(-pi,pi,largen);
    Bunit = Bunit .* exp(1i*k*freqaxis.^2);
end
  

% If necessary, scale down the B polynomial
if max(abs(Bunit)) > 1
   scaling_factor = max(abs(Bunit));
   disp(['B is scaled down by '...
       num2str(scaling_factor)]) 
   Bunit = Bunit / scaling_factor;
end

% Obtain the minimum-phase A polynomial that satisfies the 
% constraint on the norm
rhoAunit = 1 - conj(Bunit).*Bunit;
absAunit = real(sqrt(rhoAunit));
phaseAunit = findminphase(rhoAunit);
Aunit = absAunit.*exp(1i*phaseAunit);


% Check that all is well with the norm
if max(abs((conj(Bunit).*Bunit+conj(Aunit).*Aunit)-1)) > 5*eps
    display('The total norm is not 1')
end       

% Retrieve the A and B polynomial
along = ifft(fftshift(Aunit));        
blong = ifft(fftshift(Bunit));
a = along(1:polyn);
b = blong(largen/2-polyn/2+1:largen/2+polyn/2);

end