function varargout = load_targlist(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = load_target_file(varargin{:});
else
  [varargout{1}] = load_target_file(varargin{:});
end
