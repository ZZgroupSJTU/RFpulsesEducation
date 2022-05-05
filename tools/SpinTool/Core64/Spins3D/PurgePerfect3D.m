function spinsOut = PurgePerfect3D(spinsIn, spoilingDuration)
% Description: artificially sets to 0 the transverse components of the
% mangetization.

spinsOut=spinsIn;
spinsOut.M(:,:,:,1) = 0;
spinsOut.M(:,:,:,2) = 0;

if (nargin>1)
    spinsOut = Delay3D(spinsOut, spoilingDuration);
end