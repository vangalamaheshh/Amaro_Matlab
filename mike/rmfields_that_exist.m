function varargout = rmfields_that_exist(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = rmfield_if_exist(varargin{:});
else
  [varargout{1}] = rmfield_if_exist(varargin{:});
end
