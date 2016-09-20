function C = compute_mutation_rates(M,P)

% default parameter values

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'include_silent_in_sort',true);
P=impose_default_value(P,'exclude_these_genes_from_silent_calculation',{});
P=impose_default_value(P,'exclude_top_n_genes',40);
P=impose_default_value(P,'manual_correction_TSP_WashU',false);

C = [];
C.nc = length(M.mutclass)-M.NUM_INDEL_CLASSES-1;
C.nsil = sum(M.n_silent(:,1:C.nc,:),3);
C.nnon = sum(M.n_nonsilent(:,1:C.nc,:) - M.n_indel(:,1:C.nc,:),3);
C.Ntot = sum(M.N_cov(:,1:C.nc,:),3);

C.silent = M.breakdown.frac(:,1:C.nc);
C.nonsilent = M.breakdown.frac(:,C.nc+[1:C.nc]) + M.breakdown.frac(:,2*C.nc+[1:C.nc]);
C.silent_frac = C.silent ./ (C.silent + C.nonsilent);
C.nonsilent_frac = C.nonsilent ./ (C.silent + C.nonsilent);

% fix NaN's

C.silent_frac(isnan(C.silent_frac)) = 0;
C.nonsilent_frac(isnan(C.nonsilent_frac)) = 0;

C.Nsil = round(C.Ntot.*C.silent_frac);
C.Nnon = round(C.Ntot.*C.nonsilent_frac);

% append total as last column

C.Nsil = [C.Nsil sum(C.Nsil,2)];
C.Nnon = [C.Nnon sum(C.Nnon,2)];
C.nsil = [C.nsil sum(C.nsil,2)];
C.nnon = [C.nnon sum(C.nnon,2)];
C.nc = C.nc + 1;

% sort by mutations-per-MB, with most-mutated at end

if P.include_silent_in_sort
  numer = (C.nnon(:,end)+C.nsil(:,end));
  denom = (C.Nnon(:,end)+C.Nsil(:,end));
else
  numer = C.nnon(:,end);
  denom = C.Nnon(:,end);
end
denom(denom==0) = NaN;
C.gene_rate = numer ./ denom;

% save unsorted versions

C.unsorted.Nsil = C.Nsil;
C.unsorted.Nnon = C.Nnon;
C.unsorted.nsil = C.nsil;
C.unsorted.nnon = C.nnon;

% sort

[tmp j] = sort(C.gene_rate);
C.Nsil = C.Nsil(j,:);
C.Nnon = C.Nnon(j,:);
C.nsil = C.nsil(j,:);
C.nnon = C.nnon(j,:);

% add up cumulative totals
C.CNsil=zeros(M.ng,C.nc);
C.CNnon=zeros(M.ng,C.nc);
C.Cnsil=zeros(M.ng,C.nc);
C.Cnnon=zeros(M.ng,C.nc);
for g=1:M.ng
  C.CNsil(g,:) = sum(C.Nsil(1:g,:),1);
  C.CNnon(g,:) = sum(C.Nnon(1:g,:),1);
  C.Cnsil(g,:) = sum(C.nsil(1:g,:),1);
  C.Cnnon(g,:) = sum(C.nnon(1:g,:),1);
end

% compute rates
C.HAT=2;
C.LOW=1;
C.HIGH=3;
C.Rnon = zeros(C.nc,3,M.ng);
C.Rsil = zeros(C.nc,3,M.ng);

for c=1:C.nc
  [phat pci] = binofit(C.Cnnon(:,c), C.CNnon(:,c));
  C.Rnon(c,C.HAT,:) = phat;
  C.Rnon(c,C.LOW,:) = pci(:,1);
  C.Rnon(c,C.HIGH,:) = pci(:,2);
  [phat pci] = binofit(C.Cnsil(:,c), C.CNsil(:,c));
  C.Rsil(c,C.HAT,:) = phat;
  C.Rsil(c,C.LOW,:) = pci(:,1);
  C.Rsil(c,C.HIGH,:) = pci(:,2);
end

tot = repmat(C.Rnon(C.nc,C.HAT,:),[C.nc 3 1]);
z = find(~tot);
tot(z) = NaN;
C.RRnon = C.Rnon ./ tot;
C.RRnon(z) = NaN;
tot = repmat(C.Rsil(C.nc,C.HAT,:),[C.nc 3 1]);
z = find(~tot);
tot(z) = NaN;
C.RRsil = C.Rsil ./ tot;
C.RRsil(z) = NaN;

% prepare silent data

genes_for_silent = (setdiff(1:M.ng, ...
  listmap(P.exclude_these_genes_from_silent_calculation, M.gene.name(j))));

C.TNsil = sum(C.Nsil(genes_for_silent,:));
C.TNnon = sum(C.Nnon(genes_for_silent,:));
C.Tnsil = sum(C.nsil(genes_for_silent,:));

if P.manual_correction_TSP_WashU
%%%%%%%%%%%%%%%%%%%%%%2008-08-05--THIS LINE WAS CAUSING THE ERROR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  C.TNsil = round(C.TNsil * (9367499 / (C.TNsil(end)+C.TNnon(end))));   % WashU coverage totals only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% corrected version:
%  C.TNsil = round(C.TNsil * (93967499 / (C.TNsil(end)+C.TNnon(end))));   % WashU coverage totals only
%%%%
end

[phat pci] = binofit(C.Tnsil, C.TNsil);
C.RTsil = [pci(:,1) phat' pci(:,2)];
tot = repmat(C.RTsil(C.nc,C.HAT),[C.nc 3]);
C.RRTsil = C.RTsil ./ tot;

% prepare "pooled" data (silent + non-top-X nonsilent)

num_to_exclude = P.exclude_top_n_genes;

C.nontop = j(1:end-num_to_exclude);

if num_to_exclude<M.ng
  C.npool = C.Cnnon(M.ng-num_to_exclude,:);
  C.Npool = C.CNnon(M.ng-num_to_exclude,:);
else
  C.npool = 0;
  C.Npool = 0;
end

C.npool = C.npool + C.Tnsil;
C.Npool = C.Npool + C.TNsil;

[phat pci] = binofit(C.npool, C.Npool);
C.Rpool = [pci(:,1) phat' pci(:,2)];
tot = repmat(C.Rpool(C.nc,C.HAT),[C.nc 3]);
C.RRpool = C.Rpool ./ tot;

end

