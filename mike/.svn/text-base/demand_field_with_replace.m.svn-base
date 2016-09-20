function varargout = demand_field_with_replace(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = require_fields_with_convert(varargin{:});
else
  [varargout{1}] = require_fields_with_convert(varargin{:});
end
