function spins = ReturnSpinsToThermalEquilibrium(spins, initMag)
% Returns spins to thermal equilibrium along the z-axis.
%   spins = ReturnSpinsToThermalEquilibrium(spins) Works for non J-coupled
%   spins (with and without exchange)
%
%   spins = ReturnSpinsToThermalEquilibrium(spins, initMag)  Allows the
%   user to specify the initial magnetization. For non-exchanging spins,
%   this is a 3x1 vector. For exchanging spins, this can be either a 
%   3x1 vector (which will then be enforced on all chemical sites), or
%   a 3Nx1 vector, where N is the number of chemical sites. 

numSpins = numel(spins);
if nargin<2
    isInitMag = 0;
else
    isInitMag = 1;
    % Check input of initMag. Make sure it is column mode
    if isrow(initMag), initMag = initMag.'; end
    N = size(initMag, 1);
    if mod(N,3)~=0, error('initMag must have 3 elements (or a multiple of, for exchanging spins). Currently has %d elements.', N); end
end


if isfield(spins, 'K')
    % Exchanging spins with more than 1 site
    for curSpin = 1:numSpins
        numSites = numel(spins(curSpin).cs);
        
        if isInitMag
            if N==(3*numSites)
                spins(curSpin).M = initMag;
            else
                spins(curSpin).M = repmat(initMag, [numSites, 1]);
            end
        else
            % Use the equilibrium magnetization specified in M0
            X = zeros(numSites*3, 1);
            M0 = spins(curSpin).M0;
            if numel(M0)~=numSites
                M0 = M0*ones(numSites, 1);
            end
            X(3:3:end) = M0; 
            spins(curSpin).M = X;
        end
    end
else
    % Non-exchanging spins.
    if (numel(spins)==1) && (ndims(spins.M)==4)
        % A spins3D structure (using InitSpins3D)
        if isInitMag
            Nx = numel(spins.xVec);
            Ny = numel(spins.yVec);
            Nz = numel(spins.zVec);
            spins.M = repmat(shiftdim(initMag, -3),[Nx Ny Nz 1]);
        else
            spins.M(:,:,:,1) = 0;
            spins.M(:,:,:,2) = 0;
            spins.M(:,:,:,3) = spins.M0;
        end
    else
        % A "conventional" spin array structure
        for curSpin = 1:numSpins
            if isInitMag
                spins(curSpin).M = initMag;
            else
                spins(curSpin).M = [0; 0; spins(curSpin).M0];
            end
        end
    end
end