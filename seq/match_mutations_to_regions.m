function R = match_mutations_to_regions(R,X)
% R = match_mutations_to_regions(R,X)
%
% given a list of mutations (X) and a list of genomic regions (R),
%   finds which mutations are in which regions, and annotates R accoringly.
%
% NOTE: assumes both R and X are sorted by position!
%
% Mike Lawrence 2010-04-29

demand_fields(R,{'chr','start','end'});
demand_fields(X,{'indiv','chr','start','end'});

fprintf('Matching: ');
nx = slength(X);
nr = slength(R);
[iu ii ij] = unique(X.indiv);
ni = length(iu);
n = zeros(nr,ni);
r = 1;
for x=1:nx, if ~mod(x,10000), fprintf('%d/%d ',x,nx); end
  if isnan(X.chr(x)), continue; end
  while (X.chr(x)>R.chr(r) || (X.chr(x)==R.chr(r) && X.start(x)>R.end(r))) && r<nr, r=r+1; end
  if X.start(x)>=R.start(r) && X.start(x)<=R.end(r)
    n(r,ij(x)) = n(r,ij(x)) + 1;
  end
end, fprintf('\n');
R.nmuts = sum(n,2);
R.nsamps = sum(n>0,2);

fprintf('Adding details: ');
R.details = repmat({''},nr,1);
for i=1:nr, if ~mod(i,100000), fprintf('%d/%d ' ,i,nr); end
  if R.nsamps(i)==0, continue; end
  txt = '';
  z = 0;
  for j=1:ni
    if n(i,j)>0
      if ~isempty(txt), txt = [txt ',']; end
      txt = [txt iu{j} ':' num2str(n(i,j))];
    else
      z=z+1;
    end
  end
  R.details{i,1} = txt;
end,fprintf('\n');
