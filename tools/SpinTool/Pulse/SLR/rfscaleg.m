function rfs = rfscaleg(rf,t);

%  rfs = rfscaleg(rf,t)
%
%    rf  -- rf waveform, scaled so sum(rf) = flip angle
%    t   -- duration of RF pulse in ms
%    rfs -- rf waveform scaled to Gauss
%

gamma = 2*pi*4.257; % kHz/G
dt = t/length(rf);
rfs = rf/(gamma*dt);

