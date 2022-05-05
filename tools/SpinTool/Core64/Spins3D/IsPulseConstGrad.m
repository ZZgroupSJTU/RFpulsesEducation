function y = IsPulseConstGrad(pulse)
% Returns true if pulse has constant gradients (including zero gradients).

if numel(pulse.Gx) == 1
    y = true;
else
    threshold = 1e-10;
    isConstGx = sum(abs(diff(pulse.Gx)<threshold));
    isConstGy = sum(abs(diff(pulse.Gy)<threshold));
    isConstGz = sum(abs(diff(pulse.Gz)<threshold));
    if isConstGx && isConstGy && isConstGz
        y = true;
    else
        y = false;
    end
end