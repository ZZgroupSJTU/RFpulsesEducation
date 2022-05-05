function [pulseSPSP, pulseSpat, pulseSpec, outSPSP, pulseDec] = PulseCreateSPSP(SPSP)
% Syntax: [pulseSPSP, pulseSpat, pulseSpec] = PulseCreate_SPSP(SPSP)
% Creates a spatial spectral pulse according to the information stored
% in the structure SPSP, having the following fields:
%
% Field              Units     Type         Range                      Description
% amp                -         NxN double   [0,1]                      Amplitude matrix of SPSP profile
% ph                 rad       NxN double   [0,2pi]                    Phase matrix of SPSP profile
% offsets            kHz       1xN double   (-swSpec/2, swSpec/2)      Spectral offsets of peaks
% widthSpec          kHz       double       (0,swSpec)                 Width of spectral peak
% swSpec             kHz       double       (0,inf)                    Spectral width of spectral-domain
% shapeSpec          -         string        -                         Spectral pulse shape - see PulseCreate_Shaped for options
% shapeSpat          -         string        -                         Spatial pulse shape - see PulseCreate_Shaped for options
% facSpat            -         double       (0,inf)                    A "stretching" factor for spatial lobe, >~1
% FOV                mm        double       (0,inf)                    Spatial field of view (determines Ge)
% decOn              -         boolean      true, false                Whether to add a decoupling delay
% decTp              ms        double       [0,inf)                    Length of decoupling delay added (if decOn = true)
% decPwr             dB        double       (-inf, inf)                Actual power used for decoupling (with time = decTp)
% dec90              ms        double       (0, inf)                   Time for a 90-pulse on the decoupled nucleus
% decPwr90           dB        double       (-inf, inf)                Power in dB for 90-pulse on the decoupler
% decNumSteps        -         integer      1,2,3,4,...                Number of pulse steps for decoupling 180
% decGating          ms        double       (0,inf)                    Minimal time spent on gating before 180. Same time is spent after.
% decTime180         ms        double       (0,inf)                    Minimal time spent on doing the 180 decoupling pulse.
% tpSpat             ms        double       (0,inf)                    Time per "zig" - not including decoupling!
% numSpec            -         integer      1,2,3,4,...                Number of "zigzags" in our pulse
% decActualGating    ms        double       (0,inf)                    Actual time spent on gating before 180
% decActualTime180   ms        double       (0,inf)                    Actual time spent on 180 decoupling pulse
% decNumGatingSteps  -         integer      1,2,3,4,...                Number of time steps during gating
%
% NOTE: The last 4 fields are not initialized upon input. They are computed during the function's executation
% and returned to the caller.
%
% Output:
%   pulseSPSP   The full spatial-spectral pulse.
%   pulseSpat   The spatial excitation lobes (all with amp 1, phase 0).
%   pulseSpec   A structure of spectral pulses, such that pulseSpec{k} is
%               spectral polychromatic pulse corresponding to the k-th row
%               of the amplitude-phase matrix.
%   SPSP        Information for decoupled pulses
%   pulseDec    Pulse on decoupler. Note: uses the same time step as SPSP pulse.

N = length(SPSP.offsets);
pulseDec = PulseCreate_Null;

% ------------------------------------------------------------------------
%                           (A.) Create Spectral Pulses
% ------------------------------------------------------------------------

for k=1:N
    pulseSpec{k} = PulseCreate_PolyChromatic(SPSP.shapeSpec, SPSP.widthSpec, SPSP.swSpec, SPSP.offsets, SPSP.ph(k,:), SPSP.amp(k,:));
end;

tpSpec = pulseSpec{1}.tp;
nSpec = length(pulseSpec{1}.RFamp);

% ------------------------------------------------------------------------
%                       (B.) Compute Spatial Parameters
% ------------------------------------------------------------------------
    
tpSpat = tpSpec/2/nSpec;   % ms
[FWHH, FWAB, Nfac] = PulseInfo(SPSP.shapeSpat);
widthLobe = FWHH/tpSpat;
widthSpat = widthLobe * SPSP.facSpat;
swSpat = N * widthSpat;   % kHz
Ge = swSpat/SPSP.FOV;   % kHz/mm

% Compute offsets of lobes - i.e., where their centers should appear
offsetsSpat = [-swSpat/2 + widthSpat/2:widthSpat:swSpat/2 - widthSpat/2];

% Generate basic excitation pulse - used mainly for debugging
pulseSpat = PulseCreate_PolyChromatic(SPSP.shapeSpat, widthLobe, swSpat, offsetsSpat, zeros(1,N), ones(1,N));

% Calculate spatial pulse weighting 
ampSpat = zeros(N,nSpec);
phSpat  = zeros(N,nSpec);
for k=1:nSpec
    for p=1:N
        ampSpat(p,k) = pulseSpec{p}.RFamp(k); 
        phSpat(p,k)  = pulseSpec{p}.RFphase(k);
    end;
end;
% Normalize
ampSpat = ampSpat./(max(max(abs(ampSpat))));

% ------------------------------------------------------------------------
%                          (C.) Generate SPSP Pulse     
% ------------------------------------------------------------------------

pulseSPSP = PulseCreate_Null;

for k=1:nSpec
    pulsePlus          = PulseCreate_PolyChromatic(SPSP.shapeSpat, widthLobe, swSpat, offsetsSpat, phSpat(:,k)', ampSpat(:,k)');
    pulsePlus.RFamp    = pulsePlus.RFamp/nSpec*pi;
    nSpat              = length(pulsePlus.RFamp);
    pulsePlus.Gz       = Ge*ones(1,nSpat);

    pulseSPSP          = PulseConcat(pulseSPSP,pulsePlus);

    pulseMinus.tp      = tpSpat;
    pulseMinus.RFamp   = zeros(1,nSpat);
    pulseMinus.RFphase = zeros(1,nSpat);
    pulseMinus.Gx      = zeros(1,nSpat);
    pulseMinus.Gy      = zeros(1,nSpat);
    pulseMinus.Gz      = -Ge*ones(1,nSpat);
       
    pulseSPSP          = PulseConcat(pulseSPSP,pulseMinus);

    % Check if decoupling should be applied, and compute relevant quantities
    if (SPSP.decOn == true)
        % Time step for spatial-spectral pulse
        SPSP.dt = tpSpat/nSpat;
        % Compute the actual time the decoupling will need, in ms
        SPSP.decNumSteps = ceil(SPSP.decGating*2/SPSP.dt) + ceil(SPSP.decTime180/SPSP.dt);
        SPSP.decTp = SPSP.decNumSteps*SPSP.dt;
        % Compute number of steps the gating will take
        SPSP.decNumGatingSteps = ceil(SPSP.decGating/SPSP.dt);
        % Compute the actual gating time
        SPSP.decActualGating = SPSP.decNumGatingSteps*SPSP.dt;
        % Compute the actual 180-decoupling pulse time
        SPSP.decActualTime180 = ceil(SPSP.decTime180/SPSP.dt)*SPSP.dt;
        % Create decoupling delay on transmitter pulse
        pulseDecTemp.tp = SPSP.decTp;
        pulseDecTemp.RFamp = zeros(1,SPSP.decNumSteps);
        pulseDecTemp.RFphase = zeros(1,SPSP.decNumSteps);
        pulseDecTemp.Gx = zeros(1,SPSP.decNumSteps);
        pulseDecTemp.Gy = zeros(1,SPSP.decNumSteps);
        pulseDecTemp.Gz = zeros(1,SPSP.decNumSteps);
        pulseSPSP = PulseConcat(pulseSPSP,pulseDecTemp);
    end;
end;

% Make sure amplitudes are not negative. If they are, reverse the sign and add a pi
% to the phase (at the particular point in time it happens).
pulseSPSP = PulseToPositive(pulseSPSP);

SPSP.numSpec = nSpec;
SPSP.tpSpat = tpSpat;

% If decoupling has taken place, update data on output SPSP pulse
if (SPSP.decOn == true)
    % Compute power for 180 decoupling pulse.
    SPSP.decPwr = 20*log10(2*SPSP.dec90/SPSP.decTime180) + SPSP.decPwr90;
    % Compute actual gating time before 180-decoupling pulse
end;

outSPSP = SPSP;

% If decoupling is turned on, create the decoupling pulse.
if (SPSP.decOn == true)
    disp('Creating Decoupling Pulse ... ');
    for k=1:SPSP.numSpec;   % Loop over zigzags
        % Wait during transmitter irradiation
        pulseTemp = PulseCreate_Zero(tpSpat*2, nSpat*2);
        pulseDec = PulseConcat(pulseDec, pulseTemp);
        % Apply gating
        pulseTemp = PulseCreate_Zero(SPSP.decActualGating, SPSP.decNumGatingSteps);
        pulseDec = PulseConcat(pulseDec, pulseTemp);
        % Apply 180
        numStepsDec180 = ceil(SPSP.decTime180/SPSP.dt);
        pulseTemp = PulseCreate_Const(SPSP.decActualTime180, numStepsDec180, 1, 0);
        pulseDec = PulseConcat(pulseDec, pulseTemp);
        % Apply gating
        pulseTemp = PulseCreate_Zero(SPSP.decActualGating, SPSP.decNumGatingSteps);
        pulseDec = PulseConcat(pulseDec, pulseTemp);
    end;
end;


% ========================================================================
%                            Algorithm Outline
% ========================================================================
%
% (A.) Creating Spectral Pulses
% -----------------------------
% At the first step the spectral polychromatic pulses are created, using the
% function PulseCreate_PolyChromatic, having syntax
%
%   pulse = PulseCreate_PolyChromatic(shape, width, sw, offsets, phases, amps)
%
% The idea is to loop over the N rows of the amplitude & phase matrices and 
% create a spectral pulse for each. The pulses are stored in pulseSpec{k},
% for k=1..N
%
%
% (B.) Compute Spatial Parameters
% -------------------------------
% Next we compute Ge. Our basic spatial profile is of the form:
%
%       FWHH
%     ---------- * SPSP.facSpat 
%     tp_spatial
%     <--------->
%         ____        ____        ____        ____      
%        |    |      |    |      |    |      |    | 
%        |    |      |    |      |    |      |    | 
%       /      \    /      \    /      \    /      \
%     -----------------------------------------------> v
%
%     <  - - - - - - - Ge * FOV = swSpat - - - - - - >
%
%
% There are a total of N lobes, each having the same width. The total
% spectral width is given by Ge*FOV (where Ge is in kHz/mm), and also by:
%
%                                                 (    FWHH              )
%    (number of lobes) * (width of each lobe) = N*( ---------- * facSpat )
%                                                 ( tp_spatial           )
%
% where facSpat is a number >~ 1 allowing us the freedom to "shrink" or "expand"
% our spatial lobe. Hence
%
%         (number of lobes) * (width of each lobe)
%    Ge = ----------------------------------------
%                           FOV
% 
%
% Next, the weighting for the spatial pulses are computed (basically, the 
% spectral pulse's RF):
%
% Pulse (N steps)
%  _
% / \ [              Spectral RF 1                 ]
%  |  [              Spectral RF 2                 ]
%  |                         .
%  |                         . 
%  |                         . 
%  |  [              Spectral RF N                 ]
%    ----------------------------------------------> time (nSpec steps)
%          
% The matrix ampSpat contains the (normalized) magnitudes, while the
% matrix phSpat contains the phases (0..2pi). Thus, e.g., ampSpat(:,k) 
% contains the amplitudes of the k-th spatial pulse (total of nSpec).
%
%
% (C.) Generate the Spatial Spectral Pulse
% ----------------------------------------


