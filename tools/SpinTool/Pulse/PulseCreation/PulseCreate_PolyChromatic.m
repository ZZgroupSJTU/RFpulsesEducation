function pulse = PulseCreate_PolyChromatic(shape, width, sw, offsets, phases, amps)
% SYNTAX: PulseCreate_PolyChromatic(shape, width, sw, offsets, phases, amps)
%
% Creates and returns a polychromatic excitation pulse.
% Input parameters:
% Parameter     Units      Type     Range           Description
% shape         -          String                   Pulse excitation shape: 
%                                                   'gaussian4', 'gaussian5', 'sinc4'
% width         kHz        Double                   Peak width
% sw            kHz        Double   (0,inf)         Peak spectral width around 0 (-sw/2 to sw/2)
% offsets       kHz        Array    (-sw/2,sw/2)    Peak offset
% phases        rad        Array    [0..2pi]        Phases of peaks
% amp           -          Array    [0..1]          Scaling for the amplituded


N = length(offsets);

pulse = PulseCreate_Shaped(shape,width,offsets(1), sw, phases(1), amps(1));

for k=2:N
    pulse = PulseAdd(pulse, PulseCreate_Shaped(shape,width,offsets(k), sw, phases(k), amps(k)));
end;

   
