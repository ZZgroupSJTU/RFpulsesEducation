function y = isnucleus(n)

if ~ischar(n)
    y = false;
    return
end

if ismember(lower(n), {'1h', '2h', '13c', '14n', '15n', '17o', '19f', '31p'})
    y = true;
else
    y = false;
end
    