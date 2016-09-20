function M = dRanger3(sample,P)

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR2');
P = impose_default_value(P,'dRanger_stringent_realign_output_file','all.weird.pairs.nums.realigned.mat');
P = impose_default_value(P,'remove_duplicates',true);
P = impose_default_value(P,'window_size',2000);
P = impose_default_value(P,'randseed',1234);

basedir = '/xchip/tcga_scratch/lawrence';
fsuffix = [P.dRangerPreprocess_output_dir_suffix '/' P.dRanger_stringent_realign_output_file];
fprintf('Loading input data\n');
tmp = load([basedir '/' sample '/tumor_' fsuffix]);  T = tmp.X;
tmp = load([basedir '/' sample '/normal_' fsuffix]); N = tmp.X;

nt = size(T,1); nn = size(N,1);

keyboard

fprintf('Building matrix of pairs\n');
M = [ones(nt,1) zeros(nt,1) T(:,2:9);...
     zeros(nn,1) ones(nn,1) N(:,2:9)];
NTUM=1; NNORM=2;
CHR1=3; STR1=4; ST1=5; EN1=6;
CHR2=7; STR2=8; ST2=9; EN2=10;

%if P.remove_duplicates
%  fprintf('Removing duplicates\n');
%  [u ui uj] = unique(M(



nm = size(M,1);
clear T N;
[tmp ord] = sortrows(M(:,[CHR1 ST1 CHR2 ST2]));
M = M(ord,:);
M_orig = M;

M = M_orig;
r = P.window_size;
nm = size(M,1); gone = false(nm,1);
for c=1:24
    fprintf('\n  Chromosome %d:  ',c);
    cfirst = find(M(:,CHR1)==c,1);
    if c<24, clast = find(M(:,CHR1)>c,1)-1; else clast = nm; end
    nc = clast - cfirst + 1;
    ord = cheapperm(nc,randseed);
    lo = cfirst; hi = cfirst;
    for i=1:nc, idx = cfirst-1+ord(i);
      if ~mod(i-1,100000), fprintf('%d/%d ',i,nc); end
      if gone(idx), continue; end
      while M(idx,ST1)-M(lo,ST1)>r, lo=lo+1; end
      while hi<clast & M(hi+1,ST1)-M(idx,ST1)<=r, hi=hi+1; end
      neigh = lo:hi;
      clust = neigh(~gone(neigh) & M(neigh,CHR2)==M(idx,CHR2) &...
        abs(M(neigh,ST2)-M(idx,ST2))<=r &...
        M(neigh,STR1)==M(idx,STR1) & M(neigh,STR2)==M(idx,STR2));
      if length(clust)>1
        head = clust(1); rest = clust(2:end);
        M(head,[NTUM NNORM]) = sum(M(clust,[NTUM NNORM]));
        M(head,[ST1 ST2]) = min(M(clust,[ST1 ST2]));
        M(head,[EN1 EN2]) = max(M(clust,[EN1 EN2]));
        gone(rest) = true;
      end % endif
    end % next i
    fprintf('(->%d clusters)\n',sum(~gone(cfirst:clast)));
end % next c
M = M(~gone,:);








