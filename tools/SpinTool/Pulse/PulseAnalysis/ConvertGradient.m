function G = ConvertGradient(G, inputUnits, outputUnits, nucleus)
% Converts a given input gradient to a given output gradient. For
% example, from kHz/mm to mT/m, or from G/cm to mT/m. If nucleus is not
% defined, 1H is assumed.
%
% Input Variables
% Variable Name   Description
% G               The numerical value of the input gradient.
% inputUnits      Units of input gradient. Possible values:
%                 'mt/m'
%                 'khz/mm'
%                 'g/cm'
% outputUnits     Same as inputUnits, only for the output gradient value.
% nucleus         Optional. Assumed '1H' by default, but can also be set
%                 to '13C', '15N', '31P'. 

if nargin<4
    nucleus = '1h';
end

gmr = GetGyromagneticRatio(nucleus);

switch lower(inputUnits)
    case 'khz/mm'
        switch lower(outputUnits)
            case 'mt/m'
                G = G/gmr*1000;
            case 'g/cm'
                G = G/gmr*1000; % First, convert to mT/m
                G = 0.1*G; % Next, convert mT/m to G/cm
            otherwise % default kHz/mm 
                % Do nothing
        end
    case 'g/cm'
        switch lower(outputUnits)
            case 'mt/m'
                G = G*10;
            case 'g/cm'
                % Do nothing
            otherwise % default kHz/mm 
                G = G*10; % g/cm --> mT/m
                G = G*gmr*0.001; % mT/m --> kHz/mm
        end
    otherwise % mt/m by default
        switch lower(outputUnits)
            case 'mt/m'
                % Do nothing
            case 'g/cm'
                G = G*0.1; % mT/m --> G/cm
            otherwise % default kHz/mm 
                G = G*gmr*0.001;
        end
end