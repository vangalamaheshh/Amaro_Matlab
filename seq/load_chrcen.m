function varargout = load_chrcen(varargin)

if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = load_cen(varargin{:});
else
  [varargout{1}] = load_cen(varargin{:});
end


