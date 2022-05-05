function pulse = PulseCreateChirpseg(initialFreq, finalFreq, pulseDuration, numSteps, sampleLength, whichAxis)
% Description: The current function creates a pulse structure for 
% a chirped, linearly swept, homogeneous excitation pulse. 
%
% Inputs:
%
% Variable Name   Units   Description
% initialFreq     kHz     Initial excitation frequency
% finalFreq       kHz     Final excitation frequency
% pulseDuration   ms      Pulse duration
% numSteps        -       Number of excitation points
% sampleLength    mm      Sample length (used to set gradient). Optional
% whichAxis       char    Equal to 'x', 'y', 'z'. Optional
%
% Outputs:
%
% Variable Name   Units   Description
% pulse           -       OupulseDurationut pulse structure.

% Create time vector
excTimeVec = linspace(0,pulseDuration,numSteps);

% Calculate the sweep rate (often denoted R in many papers)
bandwidth1 = finalFreq - 0;
bandwidth2 = 0 - initialFreq;
bandwidth=(bandwidth1+bandwidth2);
sweepRate = bandwidth/pulseDuration;

% Set pulse parameters. 
pulse.tp = pulseDuration;
pulse.RFamp = wurst(40,numSteps); % Use a wurst-shaped RF amplitude to smooth the response of the pulse. 

% RFphase1=2*pi*(0*excTimeVec(1:round(numSteps/2)) + sweepRate*excTimeVec(1:round(numSteps/2)).^2 / 2);
% RFphase2= RFphase1(end)+2*pi*(initialFreq*(excTimeVec(round(numSteps/2+1):end)-pulseDuration/2)+...
%     (sweepRate*(excTimeVec(round(numSteps/2+1):end)-pulseDuration/2).^2)/ 2); 
% % pulse.RFphase=[RFphase1,RFphase2];
% pulse.RFphase=2*pi*(initialFreq*excTimeVec + sweepRate*excTimeVec.^2 / 2);
nseg=1000;
excTimeVecseg=reshape(excTimeVec,[],nseg);
Oseg=linspace(0,finalFreq,nseg+1);
for k=1:nseg
    w(k,:)=(-1)^(k-1)*Oseg(k)+0.5*sweepRate*excTimeVecseg(:,1).';
end

w=reshape(w.',1,[]);
for k=1:length(w);
    pulse.RFphase(k)=sum(w(1:k))*(excTimeVec(2)-excTimeVec(1));
end
% Quadratic phase profile
pulse.Gx = zeros(1,numSteps);
pulse.Gy = zeros(1,numSteps);
pulse.Gz = zeros(1,numSteps);

% Set chirp's amplitude for a 90-deg. excitation. For details, see
% appendix in Shrot & Frydman, JMR 172:179-190 (2005).
pulse.RFamp = pulse.RFamp./max(pulse.RFamp)*0.267*sqrt(abs(sweepRate));
pulse.t=excTimeVec;
% Set an appropriate field gradient
if nargin>4
    gradientStrength = bandwidth/sampleLength;   % kHz/mm
    switch whichAxis
        case 'x',
            pulse.Gx = gradientStrength.*ones(1,numSteps);
        case 'y',
            pulse.Gy = gradientStrength.*ones(1,numSteps);
        case 'z',
            pulse.Gz = gradientStrength.*ones(1,numSteps);
        otherwise,
            disp('Error in PulseCreateChirp: axis undefined. Aborting!');
            beep
            return
    end;
end

rfstruct.head.integral=0.2;

rfstruct.head.bandwidth     =bandwidth;
rfstruct.head.inversionBw   =0;  
rfstruct.head.modulation    ='phase';
rfstruct.head.rfFraction    =0;
rfstruct.head.type          ='SLR';
rfstruct.head.version       =1.0;

rfstruct.pulseName='chirptest';
rfstruct.pts=numSteps;
rfstruct.phase = mod(pulse.RFphase*180/pi,360);
rfstruct.amp=pulse.RFamp*1023/max(pulse.RFamp);
saveRFfile(rfstruct,pwd);
