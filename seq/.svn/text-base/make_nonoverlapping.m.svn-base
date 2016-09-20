function coords =  make_nonoverlapping(coords)
% make_nonoverlapping(coords)
%
% coords should be [st1 en1; st2 en2; st3 en3; ...]
%
% Mike Lawrence 2009-05-12

while(1)
  done = true; c2 = unique(coords,'rows');
  nc = size(coords,1);
  nc2 = size(c2,1);
  if nc2 ~= nc, coords = c2; done = false; nc = size(coords,1); end
  for i=1:nc-1, for j=i+1:nc
      if coords(i,1)<=coords(j,2) & coords(i,2)>=coords(j,1)
        coords(i,1) = min(coords(i,1),coords(j,1));
        coords(i,2) = max(coords(i,2),coords(j,2));
        coords(j,:) = [nan nan]; done = false;
  end,end,end
if done, break, end
end

coords = coords(~isnan(coords(:,1)),:);

