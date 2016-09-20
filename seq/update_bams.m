function varargout = update_bams(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = update_bamlist(varargin{:});
else
  [varargout{1}] = update_bamlist(varargin{:});
end
