function aidx = match_genelists(l1,l2)
% aidx = match_genelists(l1,l2)
% 
% allows "///" as separator

l1 = upper(l1);
l1 = regexprep(l1,'\s','');
l2 = upper(l2);
l2 = regexprep(l2,'\s','');

aidx = listmap(l1,l2);
idx = grep('///',l2,1);
for i=1:length(idx)
  q = l2{idx(i)};
  while ~isempty(q)
    j = strfind(q,'///');
    if isempty(j)
      g=q; q=[];
    else
      g=q(1:j-1); q=q(j+3:end);
    end
    g(g==' ')=[];
    ii = find(strcmp(g,l1));
    aidx(ii) = idx(i);
  end
end
