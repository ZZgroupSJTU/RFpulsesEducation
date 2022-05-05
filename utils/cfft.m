function res = cfft(x,dim)
% centered Fourier transforms for MRI: image space ---> k space
% Zhiyong Zhang (zhiyongxmu@gmail.com)

if nargin < 2
    dim = 1;
end

res = 1/sqrt(size(x,dim))*fftshift(fft(ifftshift(x,dim),[],dim),dim);