function phantom = LoadPhantomFromFile(filename)
% Loads a particular phantom from a .pha file.
% # indicates comments
% % indicates a MATLAB command line (e.g. to define a variable)


% =======================================================================
% Read data from file
% =======================================================================

% Load the file as a text file
[fid, msg] = fopen(filename, 'r');
numLines = 0;
while (~feof(fid))
    curLine = deblank(fgetl(fid));
    if (~(isempty(curLine)))
        if (curLine(1)~='#') % Lines with a # are comments: skip
            if (curLine(1)=='%') % Lines with % are variables. Evaluate
                evalin('caller',curLine(2:end));
            else % Store information in text
                numLines = numLines + 1;
                txt{numLines} = curLine;
            end
        end
    end
end
fclose(fid);


% =======================================================================
% Determines indices (in txt) of "begin" and "end" statements of structs.
% =======================================================================

beginIndices = strmatch('begin',txt);
endIndices = strmatch('end',txt);

if (length(endIndices)~=length(beginIndices))
    disp('Error in .pha file: begin/end mismatch.');
    beep
    return
end

numStructures = numel(beginIndices);


% =======================================================================
% Populate phantom
% =======================================================================


for curStruct = 1:numStructures
    curLine = beginIndices(curStruct) + 1;
    
    % Initialize default values
    phantom{curStruct}.type='cube';
    phantom{curStruct}.offset = [0 0 0];
    phantom{curStruct}.size = [1 1 1];
    phantom{curStruct}.numSpins = [1 1 1];
    phantom{curStruct}.chemShift = 0;
    phantom{curStruct}.T1 = 1000;
    phantom{curStruct}.T2 = 100;
    phantom{curStruct}.M0 = 1;
    phantom{curStruct}.eulerAngles = [0 0 0];

    % See if the loaded text has different values for some of the fields
    while (curLine<endIndices(curStruct))
        [valueType, lineValue] = ExtractLineValue(txt{curLine});
        switch (lower(valueType))
            case 'type'
                phantom{curStruct}.type = lineValue;
            case 'offset'
                phantom{curStruct}.offset = lineValue;
            case 'size'
                phantom{curStruct}.size = lineValue;
            case 'numberofspins'
                phantom{curStruct}.numSpins = lineValue;
            case 'chemicalshift'
                phantom{curStruct}.chemShift = lineValue;
            case 't1'
                phantom{curStruct}.T1 = lineValue;
            case 't2'
                phantom{curStruct}.T2 = lineValue;
            case 'm0'
                phantom{curStruct}.M0 = lineValue;
            case 'eulerangles'
                phantom{curStruct}.eulerAngles = lineValue;
            otherwise
                disp('Error: unrecognized value type in LoadPhantomFromFile.');
                beep;
                return
        end
        curLine = curLine + 1;
    end
end



% -----------------------------------------------------------------------
% Helper function: ExtractLineValue
% Given a line of the form [string] = [something], this returns the
% 'something' as the output. The 'something' can be one of 2 things:
% A set of numbers, which will be returned as an array, OR
% A string, which will be returned as a string.
% valueType holds the [string] part.
% -----------------------------------------------------------------------

function [valueType, lineValue]=ExtractLineValue(txt)
% Break string [1]=[2] into {[1], [2]}
brokenStr = deblank(strread(txt,'%s','delimiter','='));
if (numel(brokenStr)~=2)
    disp('Error: malformed .pha file. Entries must be of the form [property]=[values].');
    beep;
    return
end
valueType = brokenStr{1};
% Retain just the [2]
valueStr = brokenStr{2};
% Three scenarios are possible:
% 1. The 'value' is a number, in which case just a number is returned.
% 2. The 'value' is a set of numbers, separated by commas.
% 3. The 'value' is a string.
% These three scenarios will now be accounted for.
brokenStr = strread(valueStr,'%s','delimiter',',');

% The value is an array
if (numel(brokenStr)>1)
    lineValue = [];
    for idx=1:length(brokenStr)
        curValue = str2double(brokenStr{idx});
        % For each value in the list, add to MATLAB array 'lineValue'
        if (isnan(curValue))
            if (brokenStr{idx}(1) == '!')
                curValue = evalin('base', brokenStr{idx}(2:end));
            else
                disp('Error: malformed .pha file. Numeric entry seems to contain string.');
                beep;
                return
            end
        end
        lineValue = [lineValue, curValue];
    end
    return
end

% So it's not an array. Let's see if it's a string or a number.
curValue = str2double(brokenStr{1});
if (isnan(curValue))
    % Value is a string. 
    % If the first character of the string is "!", evaluate as a 
    % MATLAB command line expression
    if (brokenStr{1}(1)=='!')
        lineValue = evalin('base', brokenStr{1}(2:end));
    else
        lineValue = brokenStr{1};
    end
else % if curValue is a number:
    lineValue = curValue;
end

