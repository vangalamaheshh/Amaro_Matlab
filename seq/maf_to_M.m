function varargout = maf_to_M(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = maftoM(varargin{:});
else
  [varargout{1}] = maftoM(varargin{:});
end
