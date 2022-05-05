% function [aca] = b2a(bc);
%
% This function takes a b polynomial, and returns the minimum phase, 
%   minimum power a polynomial
%
% Inputs:
%   bc - beta polynomial coefficients
%
% Outputs:
%   aca - minimum phase alpha polynomial
%

function [aca] = b2a(bc);

fourierFactor = 8;

% Example: bc = [1 1 1]
n = length(bc); % n =3
% calculate minimum phase alpha
bcp = bc;  % bcp = [1 1 1]
bl = length(bc);  % bl = 3
blp = bl*fourierFactor;  % blp = 24
bcp(blp) = 0;  % bcp(24) = 0
bf = fft(bcp);  % bf = fft([1 1 1 0 0 0 0 0 0 0 ... 0])  (24 elements)

% Check maximum is not larger than 1
bfmax = max(abs(bf))
if bfmax>=1.0,                % PM can result in abs(beta)>1, not physical
  bf = bf/(1e-8 + bfmax);  %   we scale it so that abs(beta)<1 so that
end;                          %   alpha will be analytic

% afa = mag2mp(sqrt(1-bf.*conj(bf))); 
afa = sqrt(1-bf.*conj(bf));
afa = afa.*exp(1i*imag(hilbert(log(afa))));
aca = fft(afa)/blp;
aca = aca(n:-1:1);


