function varargout = calc_fdr_values(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = calc_fdr_value(varargin{:});
else
  [varargout{1}] = calc_fdr_value(varargin{:});
end
