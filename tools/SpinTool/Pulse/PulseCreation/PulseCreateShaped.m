function pulse = PulseCreateShaped(shape, bandwidth, offset, spectralWidth, pulsePhase, ampScaling)
% Description: Creates a shaped excitation pulse for a single band.
%
% Parameter     Units      Type     Range           Description
% bandwidth     kHz        Double   (0,inf)         Peak bandwidth
% offset        kHz        Double   any             Peak offset
% spectralWidth kHz        Double   (0,inf)         Peak spectral width around 0 (-spectralWidth/2 to spectralWidth/2)
% ampScaling    -          Double   [0..1]          Scaling for the amplitude
% pulsePhase    rad        Double   [0..2pi]        Pulse constant phase 
% shape         -          String                   Pulse excitation shape: 
%                                                   'gaussian4', 'gaussian5', 'sinc4'


% Make sure offset is in range (-spectralWidth/2,spectralWidth/2)
offset = mod(offset + spectralWidth/2, spectralWidth);

% Retrieve full width at half height, at base, and stepFactor (see below).
[FWHH,~,stepFactor] = PulseInfo(shape);   
pulseDuration = FWHH/bandwidth;  % ms
spectralWidth = spectralWidth / 2 + bandwidth;  % width is added to spectralWidth to ensure a peak at spectralWidth/2 won't "wrap around" to -spectralWidth/2
numSteps  = ceil(pulseDuration*spectralWidth*stepFactor);
dwellTime = pulseDuration/N;

switch lower(shape)
    case 'gaussian4',
        timeAxis = linspace(-1/2,1/2,numSteps+1);
        pulseAmplitude = exp(-(4*timeAxis(1:end-1)).^2);
    case 'gaussian5',
        timeAxis = linspace(-1/2,1/2,numSteps+1);
        pulseAmplitude = exp(-(5*timeAxis(1:end-1)).^2);
    case 'sinc4',
        timeAxis = linspace(-2,2,numSteps+1);
        pulseAmplitude = sinc(2*timeAxis(1:end-1));  % Take a look at how MATLAB defines SINC
    otherwise
        disp('Error: unrecognized shaped pulse name - aborting!');
        return
end

% Calibrate RF amp
maxB1 = 0.25*numSteps/(sum(pulseAmplitude)*pulseDuration)*1/sinc(offset*dwellTime)*ampScaling;

% Offset calibration. Term (- 2*pi*pulseDuration*center/2) is for off 0 correction. See (E.)
% Note the added pi to make the pulse real and positive at pulsePhase = 0.
timeAxis = [0:dwellTime:pulseDuration-dwellTime];

% Create pulse
pulse.tp = pulseDuration;
pulse.RFamp = pulseAmplitude.*maxB1;
pulse.RFphase = 2*pi*offset*timeAxis + pulsePhase - 2*pi*pulseDuration*offset/2 + pi;
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);




% =========================================================================
%                             Theory Outline
% =========================================================================
%
%
% (A.) Computing Tp, pulse duration
% ---------------------------------
% Each pulse shape is characterized by a Full Width at Half Height (FWHH)
% and a Full Width at Base (FWAB), defined by:
%
%     1
%     -- * FWHH = bandwidth of spectral shape
%     Tp
%
% This is used to compute Tp. Values for pulse shapes used are:
%   
%   Pulse Name     Function                x Evaluated on:    FWHH    FWAB
%   Gaussian4      exp(-(4*x)^2)           [-1/2, 1/2]        2.52    5.74
%   Gaussian5      exp(-(5*x)^2)           [-1/2, 1/2]        3.14    6.79
%   Sinc4          sinc(2*pi*x)/(2*pi*x)   [-2, 2]            7.9     11.32 
%
% (B.) Computing N, Number of Steps
% ---------------------------------
% We next determine the number of pulse steps. Note that the discrete time 
% steps used lead to a sinc modulation of our excitation profile, of the form:
%
%     modulation(w) = sinc(dt*omega/2)                       (Eq. 1)
%
% Sketch of Proof: let B1(t) be our time dependent field. Let Bi = B1(ti) 
% our sampled field. Then we can write our "sampled" field as 
%   B1(t) = sum_{k=0}^N Bk * Rect[tk, dt]
% where dt is our time step and tk is the time at which the k-th 
% Using the small tip angle approximation, compute Mxy = FT[B1(t)]. 
%
% Naively, we would choose 1/dt = sw/2 around 0.
%
%          <-- 1/dt -->
%          _          _           _
%         / \        / \         / \
%        /   \      /   \       /   \
%       /     \    /     \     /     \
%  --------+----------+-----------+------> v
%         -sw         0           sw
%
% But due to the modulation one cannot choose dt = 1/sw, since that would 
% lead a strong modulation (a peak at v=sw would not get excited at all!). 
% To resolve this, we choose
%
%   dt = 1/(stepFactor*sw/2)
%   N  = tp/dt = (tp*sw/2)*stepFactor
%
% where stepFactor is a number ~ 3, determined experimentally by trial and error
% (don't worry, I've already done the hard work for you ... ).
%
%   Pulse Name     stepFactor
%   Gaussian4      2.5
%   Gaussian5
%   Sinc4
%
% (C.) SetimeAxising the offset
% -----------------------
% This is done by setimeAxising the RF phase to 2*pi*offset*(time)
%
% (D.) Computing RF Calibration
% -----------------------------
% For the spin on resonance,
%
%   Tip angle = gamma*integral(B1(t)*dt)
%
% Hence, setimeAxising the tip angle to be pi/2, the maximal amplitude of B1 
% should be
%
%                       0.25*N
%   maxB1(kHz) = --------------------
%                integral(Bi(kHz))*tp
%
% To this we must add a correction due to the modulation in (Eq. 1), which (thanks
% to the way we chose dt to be small enough, due to stepFactor) is ~ 1. More precisely,
% the correct amplitude factor will be:
%
%                         maxB1(kHz)
%   maxB1(kHz)  --->  ------------------
%                     sinc(pi*offset*dt)  
%
% (E.) Offset Correction
% ----------------------
% The excitation profile will, in the absence of any specially added phase, accrue 
% a certain phase dependent on its offset, equal to:
%
%    2*pi*tp*offset/2      
%
% I will not offer any proof for this here. This is taken into account when 
% generating the pulse - see how the variable pulse_phase is defined.
% NOTE: even with this correction in place, our excitation profile in v will have
% a linear phase equal to 
%
%    pi*(tp-dt)*v + pi
%
% which will need to be taken into account in the post-processing phase.