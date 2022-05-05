function G=CalculateGradient(BW, VOI)
% G=CalculateGradient(BW, VOI)
% 
% Calculates the gradient G in mT/m needed to excite a given  VOI (in mm)
% for a given pulse bandwidth, BW (in kHz).

G = BW/VOI; % kHz/mm
G = G*1000/42.57; % mT/m
