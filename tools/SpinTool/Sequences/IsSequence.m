function isSeq = IsSequence(seq)
% isSeq = IsSequence(seq)
%
% Tests to see if seq is a valid sequence.

isSeq = true;

% Sequence must be a cell
if ~iscell(seq), isSeq = false; return; end

for idx=1:numel(seq)
    seqElement = seq{idx};
    if ~IsSequenceElement(seqElement), isSeq = false; return; end
end

end

function isElement = IsSequenceElement(seqElement)
    isElement = 1;
    if ispulse(seqElement)
        return
    else
        if ~iscell(seqElement)
            % If the element is not a pulse and not a cell array, it's not
            % a sequence element
            isElement = 0;
        else
            cmd = seqElement{1};
            numVars = numel(seqElement);
            switch lower(cmd)
                case 'pulse'
                    if numVars<2 || numVars>3, isElement = 0; return; end
                    if ~ispulse(seqElement{2}), isElement = 0; return; end
                    if numVars==3 % Phase cycle provided
                        if ~ismatrix(seqElement{3}), isElement = 0; return; end
                    end
                case 'delay'
                    if (numVars~=3)
                        isElement = 0;
                        return
                    end
                case 'hard'
                    if numVars~=3, isElement = 0; return; end
                case 'rect'
                    if ((numVars>5) || (numVars<4))
                        isElement = 0;
                        return
                    end
                case 'killmxy'
                    if numVars~=1, isElement = 0; return; end
                case 'killmz'
                    if numVars~=1, isElement = 0; return; end
                case 'thermal'
                    if numVars~=1, isElement = 0; return; end
                case 'purge'
                    if numVars~=5, isElement = 0; return; end
                case 'purgemoment'
                    if numVars~=5, isElement = 0; return; end
                case 'acquire'
                    switch numVars
                        case 2
                            if ~ispulse(seqElement{2}) 
                                isElement = 0; 
                                return; 
                            end
                        case 6
                            % Do nothing. It's an element.
                        otherwise
                            isElement = 0;
                            return;
                    end
                otherwise
                    isElement = 0;
                    return
            end                    
        end
    end
end

