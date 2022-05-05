function res = cifftn(x,dim)
% multiple dimentional inverse centered Fourier transforms for MRI:  k space --> image space
% Example res=cifftn(x,[1 4 5]);
% Zhiyong Zhang (zhiyongxmu@gmail.com)

if nargin < 2
    dim = 1;
end

res=x;
for m=1:length(dim)
    res = sqrt(size(res,dim(m)))*fftshift(ifft(ifftshift(res,dim(m)),[],dim(m)),dim(m));
end