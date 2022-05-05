function [pulse, excitationBandwidth] = PulseCreateSLR90(FOV, peakB1, minSlice, pulseType, isExport, sliceThickness, pulseAxis)
% SYNTAX: 
%
%    [pulse, excitationBandwidth] = PulseCreateSLR90(FOV, peakB1, ...
%                                           minSlice, pulseType, isExport, pulseAxis)
%
% Generates a linear-phase 90-degrees SLR-optimized pulse, using the
% Fourier coefficients which appear in deGraaf's book, "In Vivo NMR
% Spectroscopy", on p. 252).
%
% Input Parameters
% Name                 Units    Description
% FOV                  mm       Size of object (determines step size)
% peakB1               kHz      Maximal B1 amplitude
% minSlice             mm       Minimal excitable slice width.
% pulseType            -        One of 6 options:
%                               'linearR6'     - AM. column 1 of table
%                               'linearR12'    - AM. column 2 of table
%                               'linearR18'    - AM. column 3 of table
%                               'refocusedR6'  - FM. columns 4 & 5 of table
%                               'refocusedR12' - FM. columns 6 & 7 of table
%                               'refocusedR18' - FM. columns 8 & 9 of table
%                               The R-value is defined via BW*Duration
% isExport             0, 1     Save to .PTA file (having name 'SLR90.PTA')?
% sliceThickness       mm       Either set a slice thickness or set to []
%                               to create a frequency selective pulse
% pulseAxis            -        Either 'x', 'y',' 'z' or [] (no axis).
%
% Example: Create a 90-deg. excitation SLR pulse with peak B1 of 1 kHz:
% pulse90 = PulseCreateSLR90(100, 1.0, 1, 'linearR18)

% ========================================================================
% SLR Fourier Coefficients
% ========================================================================
%
%             SLR 90 - Linear                                    SLR 90 - Refocused
%  R=6           R=12         R=18              R=6                     R=12                   R=18
%  An            An           An          An          Bn           An          Bn          An          Bn
%  6.08         12.02        17.99             6.86                     12.80                 18.74   [1]

if nargin<1, FOV = 200; end
if nargin<2, peakB1 = 1; end
if nargin<3, minSlice = 1; end
if nargin<4, pulseType = 'linearr6'; end
if nargin<5, isExport = 0; end
if nargin<6, sliceThickness = []; end
if nargin<7, pulseAxis = []; end

SLRFourierCoeff = [...    
 0.16248152   0.07926787   0.05078526  0.15016398  0.00000000  0.08315260  0.00000000  0.05736874  0.00000000;
-0.33473283  -0.15568177  -0.10383633  0.19305354 -0.23411937  0.14400077 -0.08421706  0.10733739 -0.04101031;
 0.31861783   0.15832894   0.10238266 -0.11611928 -0.29165608  0.07895585 -0.14875552  0.08541360 -0.07763002;
-0.15415863  -0.16016458  -0.10380703 -0.20334391  0.11054030 -0.02087969 -0.16994331  0.04994930 -0.10503248;
 0.02225950   0.16450203   0.10432680  0.01999193  0.02688607 -0.13129505 -0.11648103  0.00316519 -0.11754264;
-0.00511679  -0.16193649  -0.10506775  0.00407566  0.00352382 -0.17436244  0.04634288 -0.05017866 -0.10825003;
 0.00156740   0.08371448   0.10762871 -0.00008750 -0.00249492  0.03958135  0.11688192 -0.09994613 -0.06913492;
-0.00088433  -0.01839556  -0.11074206 -0.00084618 -0.00367055  0.00929225 -0.00315922 -0.12402131  0.00610658;
 0.00013945   0.00737875   0.11078835 -0.00056449 -0.00233207  0.00437258  0.00233186 -0.07054397  0.10404452;
-0.00008296  -0.00405617  -0.06001204  0.00022508 -0.00187182  0.00231749  0.00032413  0.06727440  0.04444313;
 0.00005558   0.00266805   0.01607527  0.00048078 -0.00175806  0.00020077 -0.00127448  0.00210771 -0.00176137;
 0.00006721  -0.00158901  -0.00774033  0.00056717 -0.00156218 -0.00100345 -0.00162391  0.00331610  0.00242586;
 0.00010308   0.00105189   0.00487509  0.00061858 -0.00136616 -0.00115877 -0.00161452  0.00246314  0.00133322;
 0.00012230  -0.00065938  -0.00348665  0.00065605 -0.00119940 -0.00117903 -0.00104817  0.00131397  0.00025175;
 0.00014492   0.00036897   0.00243179  0.00066749 -0.00105887 -0.00083545 -0.00068895  0.00053659 -0.00044629;
 0.00015219  -0.00026906  -0.00180295  0.00067343 -0.00093827 -0.00044511 -0.00065681 -0.00002341 -0.00136098;
 0.00016416   0.00008059   0.00133012  0.00067005 -0.00083584 -0.00026498 -0.00080349 -0.00058156 -0.00134637;
 0.00017416  -0.00011846  -0.00097056  0.00065389 -0.00074915 -0.00026771 -0.00087878 -0.00089172 -0.00114407;
 0.00018167  -0.00003266   0.00075516  0.00064013 -0.00067090 -0.00025640 -0.00084376 -0.00099439 -0.00086665;
 0.00018828  -0.00008225  -0.00050716  0.00062404 -0.00060449 -0.00021260 -0.00082904 -0.00092942 -0.00061096;
 0.00019865  -0.00005674   0.00042757  0.00060931 -0.00054598 -0.00018750 -0.00081233 -0.00075986 -0.00044101];
% [1] - Bandwidth (Hz) * pulse length (ms)


% BWT   - Bandwith-Time product (bandwidth taken at Mxy/M0=0.5) = R-value
% TBT   - Transition band-Time product (transition band between
%         0.05<Mxy/M0<0.095)
% B1Rel - Relative (to square pulse of equivalent duration) peak power 

switch (lower(pulseType))
    case 'linearr6'
        An = SLRFourierCoeff(:,1);
        Bn = zeros(numel(An),1);
        BWT = 6.08;
        B1Rel = 6.15;
    case 'linearr12'
        An = SLRFourierCoeff(:,2);
        Bn = zeros(numel(An),1);
        BWT = 12.02;
        B1Rel = 12.8;
    case 'linearr18'
        An = SLRFourierCoeff(:,3);
        Bn = zeros(numel(An),1);
        B1Rel = 19.7;
        BWT = 17.99;
    case 'refocusedr6'
        An = SLRFourierCoeff(:,4);
        Bn = SLRFourierCoeff(:,5);
        BWT = 6.86;
        B1Rel = 6.66;
    case 'refocusedr12'
        An = SLRFourierCoeff(:,6);
        Bn = SLRFourierCoeff(:,7);
        BWT = 12.80;
        B1Rel = 12.0;
    case 'refocusedr18'
        An = SLRFourierCoeff(:,8);
        Bn = SLRFourierCoeff(:,9);
        BWT = 18.74;
        B1Rel = 17.4;
    otherwise
        disp('Error in PulseCreateSLR90: pulseType unrecognized. Returning 1ms long zero pulse.');
        pulse = PulseCreateZero(1, 1);
end

% For a 90 square peak, 2*pi*B*T = pi/2, or B*T=1/4, or B=1/(4T)
peakB1ForSquarePulse = peakB1/B1Rel; % kHz
pulseDuration = 1/(4*(peakB1ForSquarePulse)); % ms
excitationBandwidth = BWT/pulseDuration; % kHz

% The excitation bandwidth is assumed to be no smaller than a cm.
% Thus, given the FOV, the number of steps can be calculated.
FOVBW = (FOV/minSlice)*excitationBandwidth;
dwellTime = 1/(2*FOVBW); % ms
numSteps = ceil(pulseDuration/dwellTime);
timeAxis = linspace(0, pulseDuration, numSteps);

% Create complex RF shape
B1 = An(1).*ones(1,numSteps);
for idx=2:numel(An)
    B1 = B1 ...
       + An(idx)*cos(2*pi*(idx-1)*timeAxis/pulseDuration) ...
       + Bn(idx)*sin(2*pi*(idx-1)*timeAxis/pulseDuration);
end
B1 = B1*peakB1;

% Create pulse structure
pulse.tp = pulseDuration;

pulse.RFamp = abs(B1);
pulse.RFphase = angle(B1);
pulse.Gx = zeros(1, numSteps);
pulse.Gy = zeros(1, numSteps);
pulse.Gz = zeros(1, numSteps);

% BW = gamma*G*sliceThickness
if ~isempty(sliceThickness)
    G = excitationBandwidth/sliceThickness; % kHz/mm
    switch lower(pulseAxis)
        case 'x'
            pulse.Gx = pulse.Gx + G;
        case 'y'
            pulse.Gy = pulse.Gy + G;
        case 'z'
            pulse.Gz = pulse.Gz + G;
    end
end

if (isExport)
    ExportPulseToSiemens(pulse, 'SLR90', 'SLR90', excitationBandwidth, pi/2, ['SLR Linear 90. Peak = ', num2str(peakB1),' kHz. BW = ', num2str(excitationBandwidth),' @ ',num2str(pulse.tp),' ms. FOV >= ', num2str(FOV)], minSlice);
end