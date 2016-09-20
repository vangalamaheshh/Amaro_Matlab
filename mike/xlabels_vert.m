function varargout = xlabels_vert(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = xlabel_vert(varargin{:});
else
  [varargout{1}] = xlabel_vert(varargin{:});
end
