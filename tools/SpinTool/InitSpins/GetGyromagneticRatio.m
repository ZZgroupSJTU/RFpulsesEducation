function [gmr, spin, Q] = GetGyromagneticRatio(nucleus)
% Returns the gyromagnetic ratio of the given nucleus
% 
%   gmr = GetGyromagneticRatio Returns the gyromagnetic ratio 
%   (over 2 pi) in kHz/mT of protons: 42.577 ... kHz/mT.
% 
%   [gmr, spin] = GetGyromagneticRatio  Also returns the spin of the proton,
%   spin = 1/2.
%
%   [gmr, spin] = GetGyromagneticRatio(nucleus)  nucleus is defined via a 
%   string. Possible values (case insensitive): '1H', '13C', 31P', '14N',
%   '2H', '17O'. nucleus can also be a cell array of nuclei-strings, 
%   in which case gmr and spin will both be vectors of the same number of
%   elements.
%
%   [gmr, spin, Q] = GetGyromagneticRatio(nucleus) Returns the quadrupole
%   moment of the spin. 

if nargin<1
    nucleus = '1h';
else
    nucleus = lower(nucleus);
end

if ~iscell(nucleus), nucleus = {nucleus}; end
numNuclei = numel(nucleus);

for idx=1:numNuclei
    switch nucleus{idx}
        case '1h'
            gmr(idx) = 42.5774806;
            spin(idx) = 1/2;
            Q(idx) = 0;
        case '13c'
            gmr(idx) = 10.705;
            spin(idx) = 1/2;
            Q(idx) = 0;
        case '31p'
            gmr(idx) = 17.235;
            spin(idx) = 1/2;
            Q(idx) = 0;
        case '2h'
            gmr(idx) = 6.536;
            spin(idx) = 1;
            error('Quadrupolar moment still not defined for 2H');
        case '14n'
            gmr(idx) = 19.338/2/pi;
            spin(idx) = 1;
            error('Quadrupolar moment still not defined for 14N');
        case '15n'
            gmr(idx) = -27.12/2/pi;
            spin(idx) = 1/2;
            Q(idx) = 0;
        case '17o'
            gmr(idx) = 5.772;
            spin(idx) = 5/2;
            error('Quadrupolar moment still not defined for 17O');
        otherwise
            error('Unrecognized nucleus %s \n', nucleus);
    end
end