function [TB, BW, BS] = CalcExcPulseBandwidth(pulse, frac, simBW, numResPoints, magType)
% Assuming an excitation pulse is given which excites a certain band about
% zero frequency, this finds its BW (the range of frequencies it excites to
% within 1.0 and 1.0-frac efficiency) and transition bands (TB) (from 1=frac to
% frac). It does so by solving the Bloch equation.
% 
% Input Variables
% Variable Name     Units/Range    Description
% pulse
% frac              au/[0,1]      
% simBW             kHz/positive   Total simulation bandwidth
% numResPoints      -/integers     Number of frequency points for
%                                  simulating pulse's freq. response
% magType           -/-            String, = 'mz' OR 'mxy' OR 'dualmz'
%                                  OR 'press' (case insensitive)
%
% Output Variables
% Variable Name     Units/Range    Description
% BW                kHz            Bandwidth
% TB                kHz            Transition band
% BS                kHz            Band separation
%
% The meanings of BW, TB and PD change based on the magType. They are
% illustrated in the following work of ASCII art:
%
% magType = 'mz': BS = 0
%
%       mz
%      /|\
%       |
%     1 |______________                         _________________
%       |              \                       /
%       |               \                     / 
%       |                \                   /
%       |                 \                 /
%   min +                  \_______________/
%       |             <-TB-><-----BW------>
%       |_______________________________________________________> kHz 
%
%
% magType = 'mxy': BS = 0
%
%      mxy
%      /|\
%       |
%       |                              
%       |             <-TB-><----- BW ----->
%   max +                    ______________
%       |                   /              \
%       |                  /                \
%       |                 /                  \
%       |                /                    \
%       |               /                      \
%       |              /                        \
%     0 |_____________/                          \______________> kHz 
%
%
% magType = 'dualmz':
%
%       mz
%      /|\
%       |
%     1 |____                        ______                      ___
%       |    \                      /      \                    /
%       |     \                    /        \                  /
%       |      \                  /          \                /
%       |       \                /            \              /
%   min +        \______________/              \____________/
%       |                            <-BS-><TB-><----BW---->
%       |_______________________________________________________> kHz 
%
%
%
% magType = 'press': BS = 0
%
%       mxy
%      /|\
%       |                        <---- BW ---->
%       |
%     1 |______________          ______________          ________
%       |              \        /              \        /
%       |               \      /                \      /
%       |                \    /                  \    /
%       |                 \  /                    \  /
%   min +                  \/                      \/
%       |                   <------- BW + TB ------>
%       |_______________________________________________________> kHz 
%


BS = 0;
T1 = 1e6;
T2 = 1e6; 
freqMin = -simBW/2;
freqMax = simBW/2;
if (mod(numResPoints,2)==0)
    numResPoints = numResPoints + 1;
end
switch (lower(magType))
    case 'press'
        initMag = [1; 0; 0];
    otherwise
        initMag = [0; 0; 1];
end
[~, ~, Mz, ~, magMxy, freqAxis] = CalcPulseFreqResponseWithRelax(pulse, ...
   T1, T2, freqMin, freqMax, numResPoints, initMag, 0, 0);

switch (lower(magType))
    case 'mxy'
        % First, find the total bandwidth (TB+BW). Do this by starting in the
        % middle and "walking" towards the right, until you get back to 
        % |magMxy|=frac (equal approximately to 0)
        startIdx = round(numResPoints/2);
        while ((magMxy(startIdx)>frac) && (startIdx<numResPoints))
            startIdx = startIdx+1;
        end
        totBWIdx = startIdx;

        % Now go back left until you reach the "top" (up to a frac)
        while ((magMxy(startIdx)<(max(magMxy)-frac)) && (startIdx>1))
            startIdx = startIdx-1;
        end
        BWIdx = startIdx;

        BW = freqAxis(BWIdx)*2;
        totBW = freqAxis(totBWIdx)*2;
        TB = (totBW-BW)/2;

        BW = abs(BW);
        TB = abs(TB);
    case 'mz'
        % First, find the total bandwidth (TB+BW). Do this by starting in the
        % middle and "walking" towards the right, until you get back to 
        % Mxy=frac (equal approximately to 0)
        startIdx = round(numResPoints/2);
        while ((Mz(startIdx)<((1-frac))) && (startIdx<numResPoints))
            startIdx = startIdx+1;
        end
        totBWIdx = startIdx;

        % Now go back left until you reach the "bottom" (i.e., min(Mz)+frac)
        while ((Mz(startIdx)>(min(Mz)+frac)) && (startIdx>1))
            startIdx = startIdx-1;
        end
        BWIdx = startIdx;

        BW = freqAxis(BWIdx)*2;
        totBW = freqAxis(totBWIdx)*2;
        TB = (totBW-BW)/2;

        BW = abs(BW);
        TB = abs(TB);

    case 'dualmz'
        % Starting from the middle (Assumed to be = 1), start walking to 
        % the right until you encounter the 1-frac (this is where the
        % right lobe begins)
        startIdx = round(numResPoints/2);
        while ((Mz(startIdx)>((1-frac))) && (startIdx<numResPoints))
            startIdx = startIdx+1;
        end
        BSIdx = startIdx;
        
        % Keep going right until you reach the bottom
        while ((Mz(startIdx)>(min(Mz)+frac)) && (startIdx<numResPoints))
            startIdx = startIdx+1;
        end
        TBIdxLeft = startIdx;

        % Start doing the same thing from the rightmost point
        % Keep going right until the shape starts sloping up again
        startIdx = numResPoints;
        while ((Mz(startIdx)>(min(Mz)+frac)) && (startIdx>(TBIdxLeft+2)))
            startIdx = startIdx-1;
        end
        TBIdxRight = startIdx;
        
        BS = freqAxis(BSIdx)*2;
        TB = freqAxis(TBIdxLeft) - freqAxis(BSIdx);
        BW = freqAxis(TBIdxRight) - freqAxis(TBIdxLeft);

    case 'press'
        % Start from the leftmost position (where Mxy=1) and start walking
        % right until you start sloping down (that's the beginning of the
        % TB)
        startIdx = 1;
        magMxy(startIdx)
        while ((magMxy(startIdx)>(1-frac)) && (startIdx<numResPoints))
            startIdx = startIdx+1;
        end
        TBIdx1 = startIdx;

        % Keep going right until you hit the (1-frac) again. That'll
        % define the width of the TB
        while ((magMxy(startIdx)<(1-frac)) && (startIdx<numResPoints))
            startIdx = startIdx+1;
        end
        TBIdx2 = startIdx;

        % Now find the indices of the minima in the right and left parts
        % of the frequency response (that is, in Mxy), the distance
        % between which = BW+TB
        middleIdx = round(numResPoints/2);
        leftHalf = magMxy(1:middleIdx);
        rightHalf = magMxy(middleIdx+1:end);
        [~, minIdxLeft] = min(leftHalf);
        [~, minIdxRight] = min(rightHalf);
        minIdxRight = minIdxRight + middleIdx;

        BWplusTB = freqAxis(minIdxRight) - freqAxis(minIdxLeft);
        TB = freqAxis(TBIdx2) - freqAxis(TBIdx1);
        BW = BWplusTB - TB;
        
    otherwise
        disp('Unrecognized magType in CalcExcPulseBandwidth - aborting!');
        beep;
        return
end
                                                 

