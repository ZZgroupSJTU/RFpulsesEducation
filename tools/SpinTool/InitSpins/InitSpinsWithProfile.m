function spinsOut = InitSpinsWithProfile(chemShiftVec,sampleLength,sampleProfile)
% Description: the current function creates a 1D spins construct with 
% no J coupling. The construct has chemical shifts specified in chemShiftVec 
% in Hertz (e.g., chemShiftVec=[-450 -100 800] for 3 chemical shifts at -450 Hz,
% -100 Hz and 800 Hz). The magnetization is initialized to the values 
% specified in sampleProfile (a 1D real vector), which also dictates the
% number of spins. 
%
% Input:
%
% Name             Type                Units          Description      
% chemShiftVec     1xnumChemShifts     kHz            Vector of chemical shifts
% sampleLength     double              mm             Size of sample
% sampleProfile    1xnumSpinsPerShift  a.u. [0,1]     Initial sample profile
%                                                     (along z axis)

numSpinsPerShift = length(sampleProfile);
numChemShifts = length(chemShiftVec);

positionVecZ = linspace(-sampleLength/2, sampleLength/2, numSpinsPerShift);

counter = 0;
for curChemShift=1:numChemShifts
    for curSpin=1:numSpinsPerShift
        counter = counter + 1;
        spinsOut(counter).r  = [0; 0; positionVecZ(curSpin)];
        spinsOut(counter).M  = sampleProfile(curSpin);
        spinsOut(counter).cs = chemShiftVec(curChemShift);  % in kHz!
    end;
end;
