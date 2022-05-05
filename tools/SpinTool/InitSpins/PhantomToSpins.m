function spins = PhantomToSpins(phantom)
% SYNTAX: spins = PhantomToSpins(phantom)
%
% Converts a phantom structure to a spins structure.
% The phantom is a multi-cell object, having the following fields:
%
% type          Type of current part. String. At the moment, only 'cube' 
%               and 'sphere' are supported.
% offset        Offset of center, in mm. 1x3 vector.
% size          Size of phantom, in mm. 1x3 vector
% numSpins      Number of spins along each axis. 1x3 vector.
% chemShift     Chemical shift of spins. A number.
% T1            Longitudinal relaxation, in ms.
% T2            Transverse relaxation, in ms.
% M0            Equilibrium value of spins, in arbitrary units.
% eulerAngles	A 1x3 vector containing the Euler angles, in radians.
%
% Additional fields may be present, depending on type. Currently no
% additional fields are supported.

numStructures = length(phantom);
curSpin = 0;

% ------------------------------------------------------------------------
%
% Preallocate spin structure
%
% ------------------------------------------------------------------------

totNumSpins = 0;

for curStruct = 1:numStructures
    phantomType = phantom{curStruct}.type;
    numSpins    = phantom{curStruct}.numSpins;
    size        = phantom{curStruct}.size;
    switch (phantomType)
        case 'cube'
            totNumSpins = totNumSpins + numSpins(1)*numSpins(2)*numSpins(3);
        case 'sphere'
            % Create position vectors
            % AmirS 2017.01.30 - set position to zero, if only single point
            % in given dimension. (As done later below)
            % xVec = linspace(-size(1)/2, size(1)/2, numSpins(1));
            % yVec = linspace(-size(2)/2, size(2)/2, numSpins(2));
            % zVec = linspace(-size(3)/2, size(3)/2, numSpins(3));
            if (numSpins(1)>1)
                xVec = linspace(-size(1)/2, size(1)/2, numSpins(1));
            else
                xVec = 0;
            end
            if (numSpins(2)>1)
                yVec = linspace(-size(2)/2, size(2)/2, numSpins(2));
            else
                yVec = 0;
            end
            if (numSpins(3)>1)
                zVec = linspace(-size(3)/2, size(3)/2, numSpins(3));
            else
                zVec = 0;
            end

            for cx=1:numSpins(1)
                for cy=1:numSpins(2)
                    for cz=1:numSpins(3)
                        curX = xVec(cx);
                        curY = yVec(cy);
                        % AmirS 2017.01.30 - The if below is no longer
                        % needed, because taken care of above
                        % if numSpins(3)==1, curZ=0; else, curZ = zVec(cz); % end; % AmirS, commneted
                        curZ = zVec(cz); % AmirS, new
                        if (((curX/(size(1)/2))^2 + (curY/(size(2)/2))^2 + (curZ/(size(3)/2))^2) < 1 )
                            totNumSpins = totNumSpins + 1;
                        end
                    end
                end
            end
        % ------------------------------------------
        % Unrecognized structure: abort
        % ------------------------------------------
        otherwise
            beep
            disp('Error: unrecognized phantom structure.');
            return
    end
end

spins(totNumSpins).r = [0; 0; 0];
spins(totNumSpins).M = [0; 0; 1];
spins(totNumSpins).cs = 0;
spins(totNumSpins).T1 = 1000;
spins(totNumSpins).T2 = 100;
spins(totNumSpins).M0 = 1;
spins(totNumSpins).B0 = 0;
spins(totNumSpins).B1 = 1;
spins(totNumSpins).RS = 1;



% ------------------------------------------------------------------------
%
% Convert phantom structures to spin structures
%
% ------------------------------------------------------------------------

for curStruct = 1:numStructures
    % Extract variables for easier reference (cosmetics only)
    phantomType = phantom{curStruct}.type;
    offset      = phantom{curStruct}.offset;
    size        = phantom{curStruct}.size;
    numSpins    = phantom{curStruct}.numSpins;
    chemShift   = phantom{curStruct}.chemShift;
    T1          = phantom{curStruct}.T1;
    T2          = phantom{curStruct}.T2;
    M0          = phantom{curStruct}.M0;

    % Create spatial phantom axes
    if (numSpins(1)>1)
        xAxis = linspace(-size(1)/2, size(1)/2, numSpins(1)) + offset(1);
    else
        xAxis = offset(1);
    end
    if (numSpins(2)>1)
        yAxis = linspace(-size(2)/2, size(2)/2, numSpins(2)) + offset(2);
    else
        yAxis = offset(2);
    end
    if (numSpins(3)>1)
        zAxis = linspace(-size(3)/2, size(3)/2, numSpins(3)) + offset(3);
    else
        zAxis = offset(3);
    end
    
    switch lower(phantomType)
        % ------------------------------------------
        % Cube phantom structure
        % ------------------------------------------
        case 'cube'
            for idxX=1:numSpins(1)
                for idxY=1:numSpins(2)
                    for idxZ=1:numSpins(3)
                        curSpin = curSpin + 1;
                        spins(curSpin).M = [0; 0; M0];
                        spins(curSpin).r = [xAxis(idxX), yAxis(idxY), zAxis(idxZ)];
                        spins(curSpin).cs = chemShift;
                        spins(curSpin).T1 = T1;
                        spins(curSpin).T2 = T2;
                        spins(curSpin).M0 = M0;
                        spins(curSpin).B0 = 0;
                        spins(curSpin).B1 = 1;
                        spins(curSpin).RS = 1;
                    end
                end
            end
        % ------------------------------------------
        % Sphere phantom structure
        % ------------------------------------------
        case 'sphere'
            for cx=1:numSpins(1)
                for cy=1:numSpins(2)
                    for cz=1:numSpins(3)
                      
                        % AmirS 2017.01.30 - subtraced offset from xAxis
                        % and yAxis in two lines below. This way, the
                        % condition below on curX, curY, & curZ will match
                        % the one above when setting totNumSpins.
                        % (Otherwise we end with empty valued "positions"
                        % when offset is non-zero.)
                        % curX = xAxis(cx); % AmirS Commented
                        % curY = yAxis(cy); % AmirS Commented
                        curX = xAxis(cx) - offset(1); % AmirS, new
                        curY = yAxis(cy) - offset(2); % AmirS, new
                        % AmirS 2017.01.30 - Define curZ the same way as
                        % CurX and CurY. The 'if' below is no longer needed.
                        %if numSpins(3)==1, curZ=0; else, curZ = zVec(cz); end;
                        curZ = zAxis(cz) - offset(3); % AmirS, new
                        if (((curX/(size(1)/2))^2 + (curY/(size(2)/2))^2 + (curZ/(size(3)/2))^2) < 1 )
                            curSpin = curSpin + 1;
                            spins(curSpin).M = [0; 0; M0];
                            % AmirS 2017.01.30 - No need to remove the
                            % offset, it is included in xAxis, yAxis &
                            % zAxix, and we want it.
                            % spins(curSpin).r = [xAxis(cx), yAxis(cy), zAxis(cz)] - offset; % AmirS, commented
                            spins(curSpin).r = [xAxis(cx), yAxis(cy), zAxis(cz)]; % AmirS, new
                            spins(curSpin).cs = chemShift;
                            spins(curSpin).T1 = T1;
                            spins(curSpin).T2 = T2;
                            spins(curSpin).M0 = M0;
                            spins(curSpin).B0 = 0;
                            spins(curSpin).B1 = 1;
                            spins(curSpin).RS = 1;
                        end
                    end
                end
            end
        % ------------------------------------------
        % Unrecognized structure: abort
        % ------------------------------------------
        otherwise
            beep
            disp('Error: unrecognized phantom structure.');
            return
    end
end