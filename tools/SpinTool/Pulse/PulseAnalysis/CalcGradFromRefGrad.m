function G = CalcGradFromRefGrad(refGrad, sliceThickness, pulseDuration)
% Give the reference gradient, as noted in the Siemens PTA header, and 
% slice thickness (in mm) and pulse duration (in ms), this calculates the 
% gradient during the pulse in mT/m.
% (note: to convert to kHz/mm, just multiply by GetGyromagneticRatio('1h')/1000)

G = refGrad*(10/sliceThickness)*(5.12/pulseDuration);
