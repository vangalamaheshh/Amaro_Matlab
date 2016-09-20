function [mqz1 mqz2] = dRanger_find_local_MQZ_counts(X,filestem)

% load MQZ counts
fprintf('Loading MQZ counts...\n');
C = cell(24,1);
for chr=1:24
  f = [filestem 'chr' num2str(chr) '.mqz'];
  tmp = load_struct(f,'%f%f%f%f',0);
  C{chr} = rename_field(tmp,{'col1','col2','col3','col4'},...
  {'chr','start','end','count'});
end

% apply to data

nx = slength(X);
mqz1 = zeros(nx,1); mqz2 = zeros(nx,1);
for i=1:nx
  if ~mod(i,100), fprintf('%d/%d ',i,nx); end
  idx = find(C{X.chr1(i)}.start <= X.max1(i) & C{X.chr1(i)}.end >= X.min1(i));
  if ~isempty(idx), mqz1(i) = round(mean(C{X.chr1(i)}.count(idx))); end
  idx = find(C{X.chr2(i)}.start <= X.max2(i) & C{X.chr2(i)}.end >= X.min2(i));
  if ~isempty(idx), mqz2(i) = round(mean(C{X.chr2(i)}.count(idx))); end
end
