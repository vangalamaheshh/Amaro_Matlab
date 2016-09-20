function S = assertRange(S,field,default,min,max)
%assertRange provide default value or range check 
%
%   assertRange(S,field,default,min,max)
%
% S is a 1x1 struct whose members are parameter values. If field
% does not exist in S, it is created with the default value. If min
% or max are given, the value of the field is tested to make it is
% numeric and within the range. If the value fails the test, an
% 'InvalidParam' exception is thrown.

  if ~isfield(S,field)
      S=setfield(S,field,default);
  end
  
  if exist('max','var') && ~isempty(max)
    val = getfield(S,field);
    if isnumeric(val)
      if getfield(S,field) < min
        throwAsCaller(MException('GISTIC:InvalidParam',...
              'Parameter ''%s'' less than minimum of %f allowed for function ''%s.''',...
              field, min, whocalled));
      end
        
      if getfield(S,field) > max
        throwAsCaller(MException('GISTIC:InvalidParam',...
              'Parameter ''%s'' excedes maximum of %f allowed for function ''%s''.',...
              field, max, whocalled));
      end
    else
      throwAsCaller(MException('GISTIC:InvalidParam',...
            'Parameter ''%s'' must be numeric for function ''%s.''',...
            field, whocalled));
    end
  end
        

function name = whocalled
  [stack] = dbstack(2);
  name = stack(1).file;
  if name(end-1:end) == '.m'
    name = name(1:end-2);
  end
