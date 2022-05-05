function SLRPlotUnitCircleProfile(pulse, numPoints)
% Performs a forward SLR transform and then plots the value of the
% B-polynomial in the z-plane on the unit circle, using numPoints points.

[~,B] = SLRForwardTransform(pulse);

% z-Values
z = exp(1i*linspace(0,2*pi,numPoints));

% Invert
z = z(end:-1:1);

y = polyval(B, z);

figure
subplot(2,1,1);
plot(abs(y));
subplot(2,1,2);
plot(phase(y));



