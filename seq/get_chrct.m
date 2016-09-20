function varargout = get_chrcount(varargin)

if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = get_chrcount(varargin{:});
else
  [varargout{1}] = get_chrcount(varargin{:});
end
