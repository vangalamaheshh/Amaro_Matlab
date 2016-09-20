function [Z,h] = survey_categs(categdir)
% [Z,h] = survey_categs(categdir)
%
% Z = list of category names
% h = histogram of categories across entire genome

if categdir(1)=='/'
  dr = categdir;
else
  dr = ['/xchip/cga1/lawrence/db/' categdir];
end

Z = load_struct([dr '/categs.txt']);
Z = make_numeric(Z,'num');
mn = min(Z.num); mx = max(Z.num);

fprintf('Loading categ data from %s\n',dr);
for i=1:24, fprintf('chr%d ',i);
  tmp = load([dr '/chr' num2str(i) '.mat']);
  f = fieldnames(tmp);
  C = getfield(tmp,f{1});
  if i==1
    h = histc(C,mn:mx);
  else
    h=h+histc(C,mn:mx);
  end
end,fprintf('\n');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Z,h] = survey_categs('grizc65e');

bad = grep('bad',Z.name,1);
n = grep('N',Z.name,1);
badorn = union(bad,n);

tot = sum(h);
totbn = sum(h(badorn));



