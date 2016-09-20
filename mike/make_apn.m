function varargout = make_apn(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = makeapn(varargin{:});
else
  [varargout{1}] = makeapn(varargin{:});
end
