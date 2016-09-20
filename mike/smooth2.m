function z = smooth2(y,span,varargin)
% z = smooth2(y,span,...)
%
% same as smooth, but reduces edge effects by imbedding a background of zeros

y2 = [zeros(span,1);as_column(y);zeros(span,1)];
z2 = smooth(y2,span,varargin{:});
z = z2(span+1:span+length(y));



