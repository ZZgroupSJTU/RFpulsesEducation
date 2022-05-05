function res = cifft(x,dim)
% centered inverse Fourier transforms for MRI: k space -->image space 
% Zhiyong Zhang (zhiyongxmu@gmail.com)

if nargin < 2
    dim = 1;
end
res = sqrt(size(x,dim))*fftshift(ifft(ifftshift(x,dim),[],dim),dim);
