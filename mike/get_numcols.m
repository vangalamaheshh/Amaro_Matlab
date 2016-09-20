function varargout = get_numcols(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = get_num_cols(varargin{:});
else
  [varargout{1}] = get_num_cols(varargin{:});
end
