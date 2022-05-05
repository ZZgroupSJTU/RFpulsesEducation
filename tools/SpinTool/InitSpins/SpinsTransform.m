function spins = SpinsTransform(spins, R, x0, T)
% SpinsTransform  Performs spatial transformation on spin structure array
%   spins = SpinsTransform(spins, R)  Rotates a spin structure about the
%   origin according to the rotation matrix R (which must be orthogonal,
%   i.e. RR^T = I, and 3x3)
%
%       r -->  R*r + T
%
%   spins = SpinsTransform(spins, R, x0)  Allows specifying the rotation
%   point x0 around which the spins translate. If omitted, they will
%   be rotated about the origin.
%
%   spins = SpinsTransform(spins, R, x0, T)  Performs a post-rotation
%   translation, T, specified by the vector T. 

if nargin<3, x0 = [0, 0, 0]; end
if nargin<4, T = [0, 0, 0]; end

if isrow(x0), x0 = x0.'; end
if isrow(T), T = T.'; end

if ~isequal(size(R), [3 3]), error('Rotation matrix R must be 3x3'); end

% Check orthogonality. Some very minor deviations are acceptable.
E = R*R.' - eye(3);
if sum(abs(E(:))/9)>1e-8, error('Rotation matrix R must be orthogonal'); end

numSpins = numel(spins);

for idx=1:numSpins
    spins(idx).r = spins(idx).r - x0;
    spins(idx).r = R*spins(idx).r;
    spins(idx).r = spins(idx).r + x0 + T;
end

