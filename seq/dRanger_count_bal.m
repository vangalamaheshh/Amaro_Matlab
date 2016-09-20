function balfrac = dRanger_find_chains(sample,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'chain_threshold',1000);
P = impose_default_value(P,'show_top_n_paths',8);

X = load_struct(['/xchip/tcga_scratch/lawrence/' sample '/dRanger_results.txt']);
keyboard
% filtering criteria
X = reorder_struct(X,strcmp(X.filter,'0'));
X = make_numeric(X,{'normreads','tumreads'});X = reorder_struct(X,X.normreads==0 & X.tumreads>=4);
X = make_numeric(X,{'chr1','chr2','min1','max1','min2','max2','str1','str2','pos1','pos2'});
X.len1 = X.max1-X.min1+1;
X.len2 = X.max2-X.min2+1;
X = reorder_struct(X,X.len1>150 & X.len1<600 & X.len2>150 & X.len2<600);
%X.pos1 = round((X.min1+X.max1)/2);
%X.pos2 = round((X.min2+X.max2)/2);

% convert to easy table:   num   rearr   end   chr   pos    str    tumreads    normreads
nx = slength(X);
E = [(1:nx*2)' repmat((1:nx)',2,1) [ones(nx,1);2*ones(nx,1)] ...
   [X.chr1;X.chr2] [X.pos1;X.pos2] [X.str1;X.str2] ...
   [X.tumreads;X.tumreads] [X.normreads;X.normreads]];
ne = 2*nx;

% construct graph (adjacency matrix)
J = zeros(ne,ne);
for i=1:ne-1, for j=i+1:ne
  if E(i,2)~=E(j,2) && E(i,6)~=E(j,6) && ...
     (E(i,4)==E(j,4) & abs(E(i,5)-E(j,5))<=P.chain_threshold), J(i,j) = 1; J(j,i)=1; end
end,end

balfrac = mean(sum(J));
