function [x_joined,x_unjoined] = join_matching_rows(x)
% join_matching_rows(x)
%
% Mike Lawrence 2010-01-25

ncols = size(x,2);
nrows = size(x,1);

joined = false(nrows,1);
removed = false(nrows,1);

empty = find(x==-1);
x(empty)=nan;

x = sortrows(x);
for r=1:nrows-1
  if x(r,1)~=x(r+1,1), continue; end
  if x(r,2)~=x(r+1,2), continue; end
  if x(r,3)~=x(r+1,3), continue; end
  if x(r,4)~=x(r+1,4), continue; end
  % it's a match
  x(r,:) = max(x([r r+1],:));
  joined(r) = true;
  removed(r+1) = true;
  r=r+1;
end

x_joined = x(joined,:);
x_unjoined = x(~(joined|removed),:);





