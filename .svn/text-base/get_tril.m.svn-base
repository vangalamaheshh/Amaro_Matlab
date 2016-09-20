function [res,i,j]=get_tril(A,k)
if nargin==1
  k=-1;
end

if size(A,1)~=size(A,2)
  error('need a square matrix');
end

idx=tril(true(size(A,1)),k);
res=A(idx);

if nargout>1
 tmp=repmat((1:size(A,1))',1,size(A,1));
 i=tmp(idx);
 tmp=tmp';
 j=tmp(idx);
end

