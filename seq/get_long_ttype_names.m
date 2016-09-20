function varargout = get_long_ttype_names(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = get_long_tumor_type_names(varargin{:});
else
  [varargout{1}] = get_long_tumor_type_names(varargin{:});
end
