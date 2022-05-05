%  function [nk dt] = csg(k,mxg,mxs)
%
%  This routine takes a k-space trajectory and time warps it to
%  meet gradient amplitude and slew rate constraints.
%
%  Inputs:
%     k   --  k-space trajectory, scaled to cycles/cm
%     mxg --  maximum gradient, G/cm
%     mxs --  maximum slew rate, (G/cm)/ms
%
%  Outputs:
%     nk  --  new k-space trajectory meeting the constraints
%     dt  --  sample time for the new gradient
%
%  csg also reports the gradient duration required.  
%

%  Written by John Pauly, 1993
%
%  Oct 4, 2004 modified to use 'spline' in interp1, now that it works in
%      matlab 7

function [nk] = csg(k,mxg,mxs)

td = 1;
len = length(k);

% compute initial gradient, slew rate
g = [diff(k)]/(4.26*(td/len));
g = [g(1) g];
s = diff(g)/(td/len);
s = [s s(len-1)];

% Compute slew rate limited trajectory
ndts = sqrt(abs(s/mxs));
nt = cumsum(ndts)*td/len;
%nk = csplinx(nt,k,[1:len]*nt(len)/len);
nk = interp1(nt,k,[1:len]*nt(len)/len,'spline')';

% Apply the additional gradient amplitude constraint
g = [0 diff(nk)']'/(4.26*(nt(len)/len));
%ndtg = max(abs(g/mxg),1);
ndtg = max(abs(g),mxg);
nt = cumsum(ndtg)*nt(len)/(mxg*len);
%nk = csplinx(nt,nk,[1:len]*nt(len)/len);
nk = interp1(nt,nk,[1:len]*nt(len)/len,'spline')';

% report the waveform length
disp(sprintf('Gradient duration is %6.3f ms',nt(len)));

