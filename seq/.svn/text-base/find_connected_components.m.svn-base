function C = find_connected_components(X)
% C = find_connected_components(X)
%
% return value C has one number per record in X,
%   tells which component it belongs to.
% 
% Mike Lawrence 2010-04-28

require_fields(X,{'chr','starts','ends'});
nx = slength(X);

if ~isnumeric(X.chr)
  tmp = X.chr;
  X.chr = convert_chr(X.chr);
  idx = find(isnan(X.chr) | X.chr<1 | X.chr>24);
  if ~isempty(idx)
    fprintf('Ignoring genes on the following chromosomes:');
    count(tmp(idx));
  end
end

C = nan(nx,1);
maxc = 0;
for chr=1:24, fprintf('chr%d\n',chr);
  cidx = find(X.chr==chr);
  qi = repmat({zeros(0,3)},length(cidx),1);
  for i=1:length(cidx), if ~isempty(X.starts{cidx(i)})
    qi{i} = [X.starts{cidx(i)} zeros(length(X.starts{cidx(i)}),1) i*ones(length(X.starts{cidx(i)}),1);
             X.ends{cidx(i)} ones(length(X.ends{cidx(i)}),1) i*ones(length(X.ends{cidx(i)}),1)];
  end, end
  q = cat(1,qi{:}); q=sortrows(q);
  idx = find(q(1:end-1,2)~=1 | q(2:end,2)~=0);
  e = [q(idx,3) q(idx+1,3)];
  l = unionfind(e,length(cidx));
  C(cidx) = maxc+l;
  maxc = max(maxc+l);
end



