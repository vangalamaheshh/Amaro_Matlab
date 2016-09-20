function [E J path outtxt outmat] = dRanger_find_chains(X,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'chain_threshold',1500);
P = impose_default_value(P,'show_top_n_paths',8);
P = impose_default_value(P,'min_num_links_in_order_to_display',2);

outtxt=[];
outmat={};

if ischar(X)
  sample = X;
  X = load_struct(['/xchip/tcga_scratch/lawrence/' sample '/dRanger_results.txt']);
  % filtering criteria
  X = make_numeric(X,'normreads');X = reorder_struct(X,X.normreads<2);
  X = reorder_struct(X,~strcmp(X.filterW,'1'));
  X = make_numeric(X,{'fmapqzT1','fmapqzN1','fmapqzT2','fmapqzN2'});
  X = reorder_struct(X,X.fmapqzT1<0.2&X.fmapqzN1<0.2&X.fmapqzT2<0.2&X.fmapqzN2<0.2);
  X = reorder_struct(X,~strcmp(X.filterB,'1'));
  X = make_numeric(X,'tumreads');X = reorder_struct(X,X.tumreads>=4);
  X = make_numeric(X,{'chr1','chr2','min1','max1','min2','max2','str1','str2','pos1','pos2'});
  X.len1 = X.max1-X.min1+1;
  X.len2 = X.max2-X.min2+1;
  X = reorder_struct(X,X.len1>150 & X.len1<600 & X.len2>150 & X.len2<600);
else
  % presumes filtering already performed
  X = make_numeric(X,{'tumreads','normreads','chr1','chr2','min1','max1','min2','max2','str1','str2','pos1','pos2'});
end

if isfield(X,'individual')
  if length(unique(X.individual))>1
    error('multiple individuals not supported');
end,end


%X.pos1 = round((X.min1+X.max1)/2);
%X.pos2 = round((X.min2+X.max2)/2);

% convert to easy table:   num   rearr   end   chr   pos    str    tumreads    normreads
nx = slength(X);
E = [(1:nx*2)' repmat((1:nx)',2,1) [ones(nx,1);2*ones(nx,1)] ...
   [X.chr1;X.chr2] [X.pos1;X.pos2] [X.str1;X.str2] ...
   [X.tumreads;X.tumreads] [X.normreads;X.normreads]];
ne = 2*nx;

% construct graph (adjacency matrix)
J = sparse(ne,ne);
for i=1:ne-1, for j=i+1:ne
  if E(i,2)==E(j,2) || (E(i,4)==E(j,4) & abs(E(i,5)-E(j,5))<=P.chain_threshold), J(i,j) = 1; J(j,i)=1; end
end,end

% find paths
path = cell(ne,1);
for i=1:ne, path{i} = graphtraverse(J,i); end

% collapse to unique paths
ps = cell(ne,1);
for i=1:ne, ps{i} = sprintf('%d+',sort(path{i})); end
[tmp idx] = unique(ps);
path = path(idx);
len = zeros(length(path),1);
for i=1:length(path), len(i) = length(path{i}); end

% show top n paths
[tmp ord] = sort(len,'descend');
showed_nothing = true;
for i=1:P.show_top_n_paths
  if i>length(ord) || len(ord(i))<P.min_num_links_in_order_to_display, continue; end
  showed_nothing = false;
  fprintf2('Top path #%d:   %d links   (%d rearrangements)\n',i,len(ord(i)),len(ord(i))/2);
%  disp(E(path{ord(i)},:));
  tmp = E(path{ord(i)},:);
  fprintf2('%5s  %5s  %3s  %4s  %10s  %3s  %5s  %3s  %s\n',...
    'bkpt','rearr','end','chr','pos','str','#T','#N','site');
  outmat{end+1,1} = [];
  outmat{end}.bkpt = [];
  outmat{end}.rearr = [];
  outmat{end}.end = [];
  outmat{end}.chr = [];
  outmat{end}.pos = [];
  outmat{end}.str = [];
  outmat{end}.tumreads = [];
  outmat{end}.normreads = [];
  outmat{end}.site = {};
  outmat{end}.gene = {};
  for j=1:size(tmp,1)
    if tmp(j,3)==1
      site = X.site1{tmp(j,2)}; gene = X.gene1{tmp(j,2)};
    else
      site = X.site2{tmp(j,2)}; gene = X.gene2{tmp(j,2)};
    end
    fprintf2('%5d  %5d  %3d  %4d  %10d  %3d  %5d  %3d  %s\n',tmp(j,:),site);
    outmat{end}.bkpt(end+1,1) = tmp(j,1);
    outmat{end}.rearr(end+1,1) = tmp(j,2);
    outmat{end}.end(end+1,1) = tmp(j,3);
    outmat{end}.chr(end+1,1) = tmp(j,4);
    outmat{end}.pos(end+1,1) = tmp(j,5);
    outmat{end}.str(end+1,1) = tmp(j,6);
    outmat{end}.tumreads(end+1,1) = tmp(j,7);
    outmat{end}.normreads(end+1,1) = tmp(j,8);
    outmat{end}.site{end+1,1} = site;
    outmat{end}.gene{end+1,1} = gene;
  end

  is_circ = J(path{ord(i)}(1),path{ord(i)}(end));
  if is_circ, fprintf2('(circular)\n'); end
  fprintf2('\n\n');
end

if showed_nothing
  fprintf2('No paths with at least %d links\n',P.min_num_links_in_order_to_display);
end

  function fprintf2(varargin)
    t = sprintf(varargin{:});
    outtxt = [outtxt t];
    fprintf(varargin{:});
  end

end % main function
