function res = cfftn(x,dim)
% multiple dimentional centered Fourier transforms for MRI:  image space ---> k space
% Example res=cfftn(x,[1 4 5]);
% Zhiyong Zhang (zhiyongxmu@gmail.com)

if nargin < 2
    dim = 1;
end

res=x;
for m=1:length(dim)
    res = 1/sqrt(size(res,dim(m)))*fftshift(fft(ifftshift(res,dim(m)),[],dim(m)),dim(m));
end