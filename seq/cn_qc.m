function L = cn_qc(Z,L,R,Cr,P)
% cn_qc(Z,L,R,Cr,P)
%
% copynumber-based quality control
%
% inputs:
%
%            Z   segments
%             .sample                full sample name from aCGH output file
%             .sample_to_match       name to cross reference to Z
%             .istum                 (0 or 1)
%             .chr
%             .start
%             .end
%             .nprobes               number of probes in this segment
%             .segmean               log2 copy ratio from aCGH
%
%            L   lane info
%             .SM                    full sample name from Picard files
%             .sample_to_match       name to cross reference to Z
%             .PU                    flowcell+lane
%             .baitset               name of baitset
%             .istum                 (0 or 1)
%             .lane                  0-based lane index for this sample+TN
%
%            R   region info
%             .gene
%             .chr
%             .start
%             .end
%             .len
%             .membership            1=WE, 2=C6K+WE, 3=C2K+C6K+WE
%
%            Cr  coverage (regions) by lane
%             = |R|x|L| matrix of unnormalized coverage counts from RegionCovByLane
%
% outputs:
%            L   lane info with added columns:
%             .nreads                total number of reads in lane
%
% Mike Lawrence 2009-10-05

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'min_nprobes',50);
P = impose_default_value(P,'min_lane_nreads',10000);

% collapse region coverage to genes
fprintf('Collapsing region coverage to genes: ');
G = [];
[G.name gi gj] = unique(R.gene);
G.chr = R.chr(gi);
G.membership = R.membership(gi);
ng = slength(G);
nl = slength(L);
Cg = zeros(ng,nl);
for i=1:ng, if ~mod(i,1000), fprintf('%d/%d ',i,ng); end
  idx = find(gj==i);
  Cg(i,:) = sum(Cr(idx,:),1);
  G.start(i,1) = min(R.start(idx));
  G.end(i,1) = max(R.end(idx));
  G.len(i,1) = sum(R.len(idx));
end, fprintf('\n');
G.mid = (G.start+G.end)/2;
% sort genes by position
[G ord] = sort_struct(G,{'chr','start'});Cg = Cg(ord,:);

%  reformat array CN segments to genes
fprintf('Mapping aCGH segments to genes: ');
Z = reorder_struct(Z,Z.nprobes>=P.min_nprobes);
S=[];
[S.sample_to_match si sj] = unique(Z.sample_to_match);
S.sample = Z.sample(si); S.istum = Z.istum(si);
ns = slength(S);
Ca = nan(ng,ns);
[gchr gchri gchrj] = unique(G.chr);
for s=1:ns, if ~mod(s,10), fprintf('%d/%d ',s,ns); end
  sidx = find(sj==s);
  for chri=1:length(gchr)
    chr = gchr(chri);
    scidx = sidx(Z.chr(sidx)==chr);
    gcidx = find(gchrj==chri);
    for j=1:length(scidx), i=scidx(j);
      gidx = gcidx(G.start(gcidx)<=Z.end(i) & G.end(gcidx)>=Z.start(i));
      Ca(gidx,s) = Z.segmean(i);
end,end,end,fprintf('\n');

% exclude very low-coverage lanes
L.nreads = sum(Cg,1)';
badlanes = find(L.nreads<P.min_lane_nreads);

% attempt automatic identification of baitset type

% reported membership
L.membership_reported = zeros(slength(L),1);
L.membership_reported(grep('2000',L.baitset,1)) = 3;
L.membership_reported(grep('6k',L.baitset,1)) = 2;
L.membership_reported(grep('whole',L.baitset,1)) = 1;

% collapse genes by geneset membership

Cm = zeros(nm,nl);
for i=1:3
  idx = find(G.membership==i);
  Cm(i,:) = sum(Cg(idx,:),1);
end
Cm = bsxfun(@rdivide,Cm,sum(Cm,1));

c2k = find(Cm(3,:)>0.8);
c6k = setdiff(find(Cm(3,:)+Cm(2,:)>0.8),c2k);
L.membership_inferred = 1*ones(nl,1);
L.membership_inferred(c6k) = 2;
L.membership_inferred(c2k) = 3;
L.membership_inferred(badlanes) = zeros;
L.mixedup_membership = (L.membership_inferred ~= L.membership_reported & L.membership_inferred>0);

% report about baitset mixups

% xcount(L.baitset,L.membership);
% mixedup = find(L.mixedup_membership)
% notmixedup = find(L.membership==L.repmem & L.repmem>0);
% tmp = randperm(length(notmixedup))';
% rep = [mixedup; tmp(1:100)];
% [tmp lord] = sort_struct(L,{'repmem'});
% lord(ismember(lord,badlanes))=[];
% [tmp gord] = sort(G.membership);
% c = bsxfun(@rdivide,Cg,sum(Cg,1));
% lidx = lord(ismember(lord,rep));
% imagesc(min(0.0001,c(gord,lidx)));

% FOR EACH MEMBERSHIP SET (C2K,C6K,WE)

for mem=1:3
  % FILTER AND NORMALIZE DATA
  lidx = find(L.membership_inferred==mem);
  gidx = find(L.membership>=mem);
  G6 = reorder_struct(G,gidx);L6 = reorder_struct(L,lidx);Ca6 = Ca(gidx,:);Cg6 = Cg(gidx,lidx);
  % normalize all lanes to 1e6 reads
  Cg6 = bsxfun(@rdivide,1e6*Cg6,L6.nreads');
  % normalize to gene length
  Cg6 = bsxfun(@rdivide,100*Cg6,G6.len);
  gtot = sum(Cg6,2);
  % (filtering under these normalizations leads to better results, by excluding poorly covered genes)
  mincovperlen = 0.3*median(gtot);
  gidx = find(gtot>=mincovperlen);
  G6 = reorder_struct(G6,gidx);Cg6 = Cg6(gidx,:);Ca6 = Ca6(gidx,:);
  % divide by median of normals
  a = median(Cg6(:,~L6.istum),2);
  Cg6 = bsxfun(@rdivide,Cg6,a);

  % ALL RANK-CORRELATIONS (ALL LANES vs. ALL ARRAYS)
  [tmp rankA] = sort(Ca6,2); [tmp rankA] = sort(rankA,2);
  [tmp rankC] = sort(Cg6,2); [tmp rankC] = sort(rankC,2);
  rankcorr2 = (1 - dist(rankA',rankC','correlation')).^2;
  [tmp ordA] = sort_struct(S,{'istum','sample_to_match'});
  [tmp ordC] = sort_struct(L6,{'istum','sample_to_match'});

  keyboard

  %imagesc(rankcorr2(ordA,ordC));
  %idx=find(diff(L6.istum(ordC)));line([idx idx],[1 slength(S)],'linewidth',2,'color','k');
  %idx=find(diff(S.istum(ordA)));line([1 slength(L6)],[idx idx],'linewidth',2,'color','k');

  % save results in L
  
end
