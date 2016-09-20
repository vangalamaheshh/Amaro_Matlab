function varargout = parseinto(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = parse_into(varargin{:});
else
  [varargout{1}] = parse_into(varargin{:});
end
