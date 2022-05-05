function elementType = ReturnSequenceElementType(sequence, elementNum)
% SYNTAX: 
%
%    seqElement = ReturnSequenceElement(sequence, elementNum)
%
% Given an MR sequence (as described in ApplySequence.m) comprised of
% delays, pulses, etc ... this function returns the type of the specified
% elemnent number ('delay', 'hard', 'pulse', etc ... ).
% In particular, if: (1) elementNum is greater than the number of elements
% in the sequence, elementType is set to 'error'; (2) the element exists
% but its type is unrecognized, elementType = 'unknown'.

if numel(sequence)<elementNum
    elementType = 'error';
    return
end

element = sequence{elementNum};

if ispulse(element)
    elementType = 'pulse';
    return
end

% If it's not a pulse, it has to be some form of string command of the
% form element={'type', param1, param2, ... }
if ~iscell(element)
    elementType = 'unknown';
    return
end

elementStr = element{1};

if ~ischar(elementStr)
    elementType='unknown';
    return
end

%   {'delay', d}           
%   Applies a delay for d milliseconds
%
%   {'hard', ang, ph}      
%   Applies a hard pulse (0.1 microsecs) with a flip angle ang (deg) 
%   and phase ph (deg)
%
%   {'rect', ang, ph, d}   
%   Applies a rect pulse with a flip angle ang (deg) phase ph (deg) and 
%   duration d (ms)
%
%   {'purge', Gx, Gy, Gz, d}
%   Applies a delay d (ms) with gradients Gx, Gy, Gz along the x, y and
%   z axes (in mT/m)
%
%   {'purgemoment', kx, ky, kz, d}
%   Applies a delay d (ms) with gradient moments kx, ky, kz along the
%   x, y, z axes (in m^(-1))
%  
%   {'killmxy'}
%   Kills all magnetization in the xy plane
%
%   {'acquire', Nt, SW, Gx, Gy, Gz}
%   Acquires an FID with Nt points, SW spectral width (kHz) and gradients
%   Gx, Gy and Gz along the x, y, z axes (in mT/m)

possibleStrs = {'hard', 'delay', 'rect', 'purge', 'purgemoment', 'killmxy', 'acquire', 'pulse', 'killmz', 'acquirefov', 'acquirebw', 'thermal'};

if ismember(lower(elementStr), possibleStrs)
    elementType = lower(elementStr);
else
    elementType = 'unknown';
end