function varargout = rmfields_if_exist(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = rmfield_if_exist(varargin{:});
else
  [varargout{1}] = rmfield_if_exist(varargin{:});
end
