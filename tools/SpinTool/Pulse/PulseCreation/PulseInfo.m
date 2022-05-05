function [FWHH, FWAB,Nfac] = PulseInfo(shape)
% Syntax: [FWHH, FWAB] = PulseInfo(shape)
%
% Input: shape contains the pulse's shape. Current options:
%   gaussian4, guassian5, sinc4
% Retrieves (from pre-computed values) the full width at half height and
% the full width at the bottom of a particular pulse shape. Also returns
% Nfac - see PulseCreateShaped for information about this parameters
% (it is only used there).

switch lower(shape)
    case 'gaussian4',
        FWAB = 5.74;
        FWHH = 2.52;
        Nfac = 2;
    case 'gaussian5',
        FWAB = 6.79;
        FWHH = 3.14;
        Nfac = 2;
    case 'sinc4',
        FWAB = 11.32;
        FWHH = 7.9;
        Nfac = 2;  % NOT VERIFIED
    otherwise
        disp('Error: unrecognized shaped pulse name - aborting!');
        return
end
