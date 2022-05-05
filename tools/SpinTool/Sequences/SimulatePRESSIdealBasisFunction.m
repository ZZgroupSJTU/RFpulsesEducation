function [spec, freqAxis] = SimulatePRESSIdealBasisFunction(metabolite, varargin)
% Simulates a PRESS metabolite basis function
%   spec = SimulatePRESSIdealBasisFunction(metabolite) Simulates a given
%   metabolite's basis function in response to an ultrashort PRESS
%   acquisition (TE1=TE2=0.1 ms) using idealized hard pulses. metabolite
%   is a string variable which can be generated using SpinsJAddMolecule,
%   e.g. 'glu', 'naa', 'cho' and so on. See SpinsJAddMolecule for a 
%   complete list. Assumes B0=2.9T, 256 ms readout, 512 points, Gaussian
%   lineshapes with FWHM of 1 Hz.
%
%   [spec, freqAxis] = SimulatePRESSIdealBasisFunction(metabolite)
%   Also returns the frequency axis in PPM, which can then be used to
%   plot the spectrum. For example: plot(freqAxis, real(spec))
% 
%   [spec, freqAxis] = SimulatePRESSIdealBasisFunction(metabolite, 'name', value, ... )
%   Allows tailoring of output with name-value pairs:
%     Sequence-related parameters:
%     B0                B0 field, in Tesla. Default: 2.9
%     TE1              
%     TE2               Echo times of PRESS sequence. Default: 0.1 for both.
%     
%     Output-related parameters:
%     freqAxisUnits     'Hz', 'PPM' (default)
%     domain            'time', 'freq' (default)
%     lineShape         'gaussian' (default), 'lorentzian', 'voigt'
%     FWHM              Full width at half max, in Hz.
%     SW                Acquisition spectral width, in Hz. Default: 1000.
%     numAcqPts         Number of acquisition points. Default: 512.

p = inputParser;
p.addParameter('TE1', 0.1, @(x) isnumeric(x)); 
p.addParameter('TE2', 0.1, @(x) isnumeric(x));
p.addParameter('B0', 2.9, @(x) isnumeric(x));
p.addParameter('SW', 2000, @(x) isnumeric(x));
p.addParameter('numAcqPts', 512, @(x) isnumeric(x));
p.addParameter('freqAxisUnits', 'ppm', @(x) ismember(lower(x), {'ppm', 'hz'}));
p.addParameter('domain', 'freq', @(x) ismember(lower(x), {'freq', 'time'}));
p.addParameter('FWHM', 1, @(x) isnumeric(x));
p.addParameter('lineShape', 'gaussian', @(x) ismember(lower(x), {'gaussian', 'voigt', 'lorentzian'}));

p.parse(varargin{:});
inputParams = p.Results;

% ========================================================================
% Define spin system
% ========================================================================

csCenter = 4.7; % ppm
B1 = 1; % Scaling
isSecular = 0;
linewidth = 0.01; % Hz

spins = InitSpinsJ(csCenter, inputParams.B0, isSecular, linewidth, B1); 
spins = SpinsJAddMolecule(spins, metabolite);

dt = 1/inputParams.SW; % Sec
dv = inputParams.SW/inputParams.numAcqPts;
timeAxis = [0:dt:(inputParams.numAcqPts-1)*dt]*1e3; % ms
freqAxis = [-inputParams.SW/2:dv:inputParams.SW/2-dv]*1e-3; % kHz

% Calculate the inter-pulse delays, from center to center
t1 = inputParams.TE1/2; 
t2 = inputParams.TE1/2 + inputParams.TE2/2;
t3 = inputParams.TE2/2;
% Check to see delays are non-negative
if t1<0, error('Negative t1 = %.2f', t1); end
if t2<0, error('Negative t2 = %.2f', t2); end
if t3<0, error('Negative t3 = %.2f', t3); end

seq = {{'hard', 90, 270},                                  {'delay', t1}, ...
       {'hard', 180, 0},                                   {'delay', t2}, ...
       {'hard', 180, 0},                                   {'delay', t3}, ...
      };

% ========================================================================
% Apply sequence
% ========================================================================

spec = [];

[spinsOut, ~] = ApplySequenceJ(spins, seq);
if ismember(lower(inputParams.domain), {'time'})
    [~, spec]  = CreateTransitionTableJ(spinsOut, 'FWHM', inputParams.FWHM, 'isCorrectDC', true, 'LineShape', inputParams.lineShape, 'Domain', 'time', 'TimeAxis', timeAxis, 'AcqType', 'separate', 'dV', 1);
    freqAxis = timeAxis;
end
if ismember(lower(inputParams.domain), {'freq'})
    [~, spec] = CreateTransitionTableJ(spinsOut, 'FWHM', inputParams.FWHM, 'isCorrectDC', true, 'LineShape', inputParams.lineShape, 'Domain', 'freq', 'FreqAxis', freqAxis, 'AcqType', 'separate', 'dV', 1);
    switch lower(inputParams.freqAxisUnits)
        case 'ppm'
            freqAxis = freqAxis/(1e-3*inputParams.B0*GetGyromagneticRatio('1h')) + 4.7;
        case 'hz'
            % Do nothing
        otherwise
            error('You are generating a basis function in the frequency domain. Your frequency axis units can be either ppm or hz (currently %s)', inputParams.freqAxisUnits);
    end
end


spec = spec{1}.';
