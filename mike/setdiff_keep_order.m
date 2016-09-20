function varargout = setdiff_keep_order(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = setdiff_keepord(varargin{:});
else
  [varargout{1}] = setdiff_keepord(varargin{:});
end
