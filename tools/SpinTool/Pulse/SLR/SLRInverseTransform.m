function pulse = SLRInverseTransform(A,B,pulseDuration)
% Maps the two SLR polynomials, A and B, to an RF pulse with a given
% duration.
% Ap(z) = A(1) + A(2)*(1/z) + A(3)*(1/z^2) + ... + A(N)*(1/z^(N-1))
% Bp(z) = B(1) + B(2)*(1/z) + B(3)*(1/z^2) + ... + B(N)*(1/z^(N-1))

% The final elements in the pulse will be numSteps-1
numSteps = length(A);

% Initialize pulse
pulse.tp = pulseDuration;
pulse.RFamp = zeros(1, numSteps);
pulse.RFphase = zeros(1, numSteps);
pulse.Gx = zeros(1, numSteps);
pulse.Gy = zeros(1, numSteps);
pulse.Gz = zeros(1, numSteps);

dwellTime = pulseDuration/numSteps;

for k=numSteps:-1:1
%     fprintf('Inverse SLR, step %d\n', k);

    phik = 2*atan(abs(B(1)/A(1)));
    phiRF = mod(angle(-1i*B(1)/A(1)), 2*pi);
    phiRF(phiRF>2*pi-1e-5) = 0;
    pulse.RFamp(k) = phik/(2*pi*dwellTime);   % phik is in radians, phik/2pi*dt is in kHz)
    pulse.RFphase(k) = phiRF;
    Ck = cos(phik/2);
    Sk = 1i*sin(phik/2)*exp(1i*phiRF);
    Atemp =  Ck*A(1:end-1) + conj(Sk)*B(1:end-1); 
    Btemp = -Sk*A(2:end)   +       Ck*B(2:end);
    A = Atemp;
    B = Btemp;

end
