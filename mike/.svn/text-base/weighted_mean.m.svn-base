function m = weighted_mean(values,weights,dim)
% weighted_mean(values,weights,dim)
%
% Mike Lawrence 2008-06-21

% if length(values)~=length(weights), error('values and weights must be vectors of equal length'); end

if ~exist('dim','var')
  if size(values,1)==1
    dim=2;
  else
    dim=1;
  end
end

if size(values,1)==1 && size(values,2)>1 && size(weights,1)==size(values,2)
  weights = weights';
end
if size(values,2)==1 && size(values,1)>1 && size(weights,2)==size(values,1)
  weights = weights';
end

a = bsxfun(@times,weights,values);
s = sum(a,dim);
m = bsxfun(@rdivide,s,sum(weights));

