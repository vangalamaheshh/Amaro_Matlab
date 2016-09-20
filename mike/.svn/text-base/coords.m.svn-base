function c = coords(sz,i)
%
% converts matrix serial indices to full coordinates
%
% e.g. 
% X = zeros(4,5,7,2);
% X(3,1,6,1) = 8;
% X(2,5,1,2) = 8;
% find(X==8)
% --> ans = [103; 158]
% coords(size(X),ans)
% --> ans = [3 1 6 1; 2 5 1 2]
%
% Mike Lawrence 2008-11-13

ndims = length(sz);
r = cumprod(sz);
r = [fliplr(r(1:end-1)) 1];

i=i-1;

for x=1:length(i)
  for d=1:ndims
    c(x,d)=floor(i(x)/r(d));
    i(x)=i(x)-c(x,d)*r(d);
  end
end

c = fliplr(c)+1;
