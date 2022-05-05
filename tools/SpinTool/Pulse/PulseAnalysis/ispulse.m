function y = ispulse(pulse)
% Returns true if 'pulse' is a pulse object, false otherwise

y=false;

if ~isstruct(pulse)
    return
end

if ~isfield(pulse, 'tp')
    return
end
    
if ~isfield(pulse, 'RFphase')
    return
end
    
if ~isfield(pulse, 'RFamp')
    return
end
    
if ~isfield(pulse, 'Gx')
    return
end

if ~isfield(pulse, 'Gy')
    return
end

if ~isfield(pulse, 'Gz')
    return
end
    
y=true;