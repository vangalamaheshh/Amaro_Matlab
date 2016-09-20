function varargout = load_struct_noheaders(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = load_struct_noheader(varargin{:});
else
  [varargout{1}] = load_struct_noheader(varargin{:});
end
