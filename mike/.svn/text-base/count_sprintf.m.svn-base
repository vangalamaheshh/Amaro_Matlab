function T = count(a, sort_by_frequency)
%% finds the unique elements of array and counts how many of each there

if size(a,1)==1
   a = a';
end

[u ui uj] = unique(a);
nu = length(u);

% make sure _a_ and _u_ are cell arrays of strings

if ~iscell(a)
  tmpa=a;
  a = cell(length(tmpa),1);
  for i=1:length(tmpa)
    if ischar(tmpa(i))
      a(i) = {tmpa(i)};
    else
      a(i) = {num2str(tmpa(i))};
    end
  end
  tmp=u;
  u = cell(length(tmp),1);
  for i=1:length(tmp)
    if ischar(tmp(i))
      u(i) = {tmp(i)};
    else
      u(i) = {num2str(tmp(i))};
    end
  end
end

b = zeros(nu,1);
for j=1:nu
    b(j) = length(find(uj==j));
end

if exist('sort_by_frequency', 'var')
  if sort_by_frequency == 1
    [b ord] = sort(b);
    u = u(ord);
  elseif sort_by_frequency == -1
    [b ord] = sort(b, 'descend');
    u = u(ord);
  end
end
  

bb = cell(nu,1);
for j=1:nu
    bb{j} = b(j);
end

u = [u; '----TOTAL'];
bb = [bb; {length(a)}];
nu=nu+1;

L=zeros(nu,1);
for i=1:nu
 L(i)=length(u{i});
end
maxL = max(L);

T = [];

T = [T sprintf('\n')];
f = ['    %' num2str(maxL) 's: [%d]\n'];
for i=1:nu
  T = [T sprintf(f, u{i}, bb{i})];
end



