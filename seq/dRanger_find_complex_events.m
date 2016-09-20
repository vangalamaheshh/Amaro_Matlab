function p = dRanger_find_complex_events(x,adjmaxdist,mingroupsize)
% Mike Lawrence 2010-01-04

if ~exist('adjmaxdist','var'), adjmaxdist = 1e4; end
if ~exist('mingroupsize','var'), mingroupsize = 3; end

require_fields(x,{'individual','name','num','chr1','pos1','chr2','pos2'});
x = make_numeric(x,{'num','chr1','pos1','chr2','pos2'});

p = {}; p.indiv=[]; p.len=[]; p.path=[]; p.mem=[]; p.chr=[]; p.num=[]; p.name=[];
[u ui uj] = unique(x.individual);
for i=1:length(u), fprintf('%d/%d ',i,length(u));
  y = reorder_struct(x,uj==i); ny = slength(y);
  d = nan(ny,ny,4);  % distances
  d(:,:,1) = bsxfun(@minus,y.pos1,y.pos1');d(:,:,2) = bsxfun(@minus,y.pos1,y.pos2');
  d(:,:,3) = bsxfun(@minus,y.pos2,y.pos1');d(:,:,4) = bsxfun(@minus,y.pos2,y.pos2');
  a = nan(ny,ny,4);
  a(:,:,1) = bsxfun(@eq,y.chr1,y.chr1');a(:,:,2) = bsxfun(@eq,y.chr1,y.chr2');
  a(:,:,3) = bsxfun(@eq,y.chr2,y.chr1');a(:,:,4) = bsxfun(@eq,y.chr2,y.chr2');
  d(~a) = nan; d = min(abs(d),[],3);
  J = sparse(d<adjmaxdist); % adjacency matrix
  path = cell(ny,1); % paths
  for j=1:ny, path{j} = graphtraverse(J,j); end
  ps = cell(ny,1);   % collapse to unique paths
  for j=1:ny, ps{j} = sprintf('%d+',sort(path{j})); end
  ps = regexprep(ps,'\+$','');
  [tmp idx] = unique(ps);  path = path(idx);
  len = zeros(length(path),1);  % path lengths
  for j=1:length(path), len(j) = length(path{j}); end
  pm = repmat({''},length(path),1);  % path members
  for j=1:length(path), for k=1:length(path{j}), z=path{j}(k);
    pm{j} = [pm{j} sprintf('[chr%d:%d/chr%d:%d]',y.chr1(z),y.pos1(z),y.chr2(z),y.pos2(z))];
  end,end
  pc = cell(length(path),1);
  for j=1:length(path), pc{j} = sprintf('%d+',unique([y.chr1(path{j});y.chr2(path{j})])); end
  pc = regexprep(pc,'\+$','');
  pnum = cell(length(path),1);
  for j=1:length(path), pnum{j} = y.num(path{j}); end
  pname = cell(length(path),1);
  for j=1:length(path), pname{j} = y.name(path{j}); end
  idx = find(len>=mingroupsize);   % store long paths
  p.indiv = [p.indiv; repmat(u(i),length(idx),1)];
  p.len = [p.len; len(idx)]; p.path = [p.path; path(idx)];
  p.mem = [p.mem; pm(idx)]; p.chr = [p.chr; pc(idx)];
  p.num = [p.num; pnum(idx)]; p.name = [p.name; pname(idx)];
end, fprintf('\n');

p = sort_struct(p,'len',-1);

