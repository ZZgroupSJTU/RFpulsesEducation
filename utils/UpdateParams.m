function [ParamsUse] = UpdateParams(ParamsDefault, ParamsIn, WarnOn)

% [ParamsUse] = UpdateParams(ParamsDefault, ParamsIn)
%
% Creates a structure of parameters, based on a default set and an
%  overriding set (can be empty or partial). I.e. each parameters defined
%  in ParamsIn overrides the value in ParamsDefault, if non is defined in
%  ParamsIn , the default value is used.
%
% Inputs
% ------
% - ParamsDefault - Structure containing default values
% - ParamsIn - Empty variable or structure containing all or part of the
%   parametrs in ParamsDefault. If parameter is defined in ParamsIn, it
%   overrides the value in ParamsDefault.
% - WarnOn - Whether (or not) to show warnings. Defaults to 'true', i.e.
%   show warnings. (Warnings are set off 
%
% Outpust
% -------
% - ParamsUse - A structure containing all the defined values, using
%   default values if they were not overridden by the ParamsIn.
%
% Notes:
% 1) Warns if a parameter is defined in the new set (ParamsIn), but not
%    defined in the default set (ParamsDefault) ;
% 2) For now(?) - Does not go recursively into sub-structures of the
%    structure. Structures should have 1 level only.

  switch nargin
    case 2
      WarnOn = true ;
    case 3
      % Do nothing
    otherwise
      error('Too many or too few input arguments') ;
  end

  ParamsUse = struct() ; % initilize an empty structyure
  
  
  if isempty(ParamsIn)
    ParamsUse = ParamsDefault ;
    return ;
  end

  if isempty(ParamsDefault)
    ParamsUse = ParamsIn ;
    if (WarnOn)
      warning(sprintf(['ParamsDefault structure (%s) is empty, using ', ...
                       'ParamsIn (%s) as is'], ...
                      inputname(1), inputname(2))) ;
    end

    return ;
  end

  FieldNames = union(fieldnames(ParamsDefault), fieldnames(ParamsIn)) ;
  
  % Go over all fields and set value of field in ParamsUse according to
  % default value and input value
  for FieldCounter = 1:length(FieldNames) 
    CurrentField = FieldNames{FieldCounter} ;
    
    ParamsUse.(CurrentField) = SetDefaultOrNew(ParamsDefault, ...
                                               ParamsIn, ...
                                               CurrentField, WarnOn) ;
  end


return ;



function [ValueUse] = SetDefaultOrNew(ParamsDefault, ParamsNew, ...
                                      FieldName, WarnOn)

  ValueUse = [] ;
   
  if isfield(ParamsDefault, FieldName) % Field exist in Default params
    
    if isstruct(ParamsDefault.(FieldName)) && ...
       isfield(ParamsNew, FieldName)
      % ParamsDefault.(FieldName) is a structure, and ParamsNew.(FieldName)
      %  exists, so we recursively run UpdateParams again, on the
      %  substructure.
      [ValueUse] = UpdateParams(ParamsDefault.(FieldName), ...
                                ParamsNew.(FieldName), WarnOn) ;
                              
    elseif (isfield(ParamsNew, FieldName) && ...
           ~isempty(ParamsNew.(FieldName)) )
       
      % Field name is defined and non-empty in new parameters
      ValueUse = ParamsNew.(FieldName) ;
         
    else
      % Field name is not defined or is empty in new parameters, but is
      %  defined in the default set.
      ValueUse = ParamsDefault.(FieldName) ;
    end
    
  else % Field does NOT exist in Default params
    if (WarnOn)
      warning(sprintf('Field ''%s'' has no default value', FieldName)) ;
    end
      
    if (isfield(ParamsNew, FieldName))
      % Field name is defined, but may be empty in new parameters (and
      % is not defined in default parameters
      ValueUse = ParamsNew.(FieldName) ;
    else
      if (WarnOn)
          warning(sprintf(['Field ''%s'' has no default or given value,' ...
                           ' using an empty value'], FieldName)) ;
      end
      ValueUse = [] ;
    end
  end
   

return ;

