function minphase = findminphase(ro)
% FINDMINPHASE yields the minimum phase for a polynomial that has the
% squared modulus rho on the unit circle. 
% Adapted from the thesis of P. Le Roux (algorithm 3, p. 70)

% The recipe says: 
% i/ take the logarithm of the modulus
% ii/ inverse Fourier tranform, multiply by the sign function and Fourier
% transform
% iii/ take the imaginary part

nfft = length(ro);
% The function is only applicable to an even number of points
if abs(nearest(nfft/2)-nfft/2) > eps
    error('The number of elements in ro should be even');
end
    
ro(ro==0)=1e-8;
% take the logarithm of the modulus
lro = 0.5 * log(ro);

% inverse Fourier transform
lro = ifft(lro);
lro = fftshift(lro);

% multiply by the sign function
lro(1:nfft/2) = -lro(1:nfft/2);
lro(nfft/2+1) = 0;

% Fourier transform
lro = fftshift(lro);
lro = fft(lro);

if any(isnan(lro))
   
    error('Something went wrong with the minimum-phase polynomial')
    
end

% take the imaginary part
minphase = imag(lro);

end