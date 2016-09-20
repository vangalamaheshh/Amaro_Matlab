function varargout = run_MutSig_S2N(varargin)
if nargout>1
  varargout = cell(nargout,1);
  [varargout{:}] = run_MutSigS2N(varargin{:});
else
  [varargout{1}] = run_MutSigS2N(varargin{:});
end
