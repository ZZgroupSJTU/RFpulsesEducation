function spins = DelayJ(spins, T, Gx, Gy, Gz)

if nargin<5
    Gx = 0;
    Gy = 0;
    Gz = 0;
end
pulse.tp = T;
pulse.RFamp = 0;
pulse.RFphase = 0;
pulse.Gx = Gx;
pulse.Gy = Gy;
pulse.Gz = Gz;

[~, spins] = PropagateJ(spins, pulse);