function varargout = preprocess_expression(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = preprocess_BMR_covariates(varargin{:});
else
  [varargout{1}] = preprocess_BMR_covariates(varargin{:});
end
