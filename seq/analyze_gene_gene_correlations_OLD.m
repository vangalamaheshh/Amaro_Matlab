function R = analyze_gene_gene_correlations(M,P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'genes_to_consider','*required*');
P = impose_default_value(P,'restrict_to_one_mut_per_gene',true);
P = impose_default_value(P,'gene_gene_correlations_num_permutations',10000);
P = impose_default_value(P,'gene_gene_correlations_report_filename',[]);
P = impose_default_value(P,'correlation_significance_threshold',0.05);
P = impose_default_value(P,'randseed',1234);

gc = listmap(P.genes_to_consider,M.gene.name);
notfound = isnan(gc);
if any(notfound)
  fprintf('Genes not found:');disp(gc(notfound));
end
gc = gc(~notfound);
ngc = length(gc);
if ngc<2
  fprintf('Need at least 2 genes to consider.\n');
  return
end

npairs = (ngc * (ngc-1)) / 2;
R = [];
R.gene1 = cell(npairs,1);
R.gene2 = cell(npairs,1);
R.len_pc_obs = nan(npairs,1);
R.Nc = nan(npairs,1);
R.Ne = nan(npairs,1);
R.c = nan(npairs,1);
R.pc_ct = nan(npairs,1);

% compute observed statistics

T = [M.mut.gene(M.use_nonsilent) M.mut.patient(M.use_nonsilent)];
if P.restrict_to_one_mut_per_gene, T = unique(T,'rows'); end
idx_obs = cell(ngc,1); for g=1:ngc, idx_obs{g} = find(T(:,1)==gc(g)); end
ridx = 1;
for g1=1:ngc-1, for g2=g1+1:ngc
  R.gene1{ridx} = M.gene.name{gc(g1)};
  R.gene2{ridx} = M.gene.name{gc(g2)};
  m1_obs = idx_obs{g1};
  m2_obs = idx_obs{g2};
  pc_obs = intersect(T(m1_obs,2),T(m2_obs,2));
  R.len_pc_obs(ridx) = length(pc_obs);
  R.Nc(ridx) = sum(ismember(T([m1_obs;m2_obs],2),pc_obs));
  R.Ne(ridx) = length([m1_obs;m2_obs])-R.Nc(ridx);
  R.c(ridx) = 0;
  R.pc_ct(ridx) = 0;
  ridx=ridx+1;
end,end

% do permutations

fprintf('Performing permutations:\n');
rand('twister',P.randseed);
nruns = P.gene_gene_correlations_num_permutations;
idx_perm = cell(ngc,1); idx_perm_length = nan(ngc,1);
tt=tic;
for i=1:nruns, if ~mod(i,100), fprintf('%d\t%d/%d\n',toc(tt),i,nruns); end
  Tp = [T(randperm(size(T,1)),1) T(:,2)];
  for g=1:ngc
    thisgene = (Tp(:,1)==gc(g));
    idx_perm{g} = find(thisgene);
    idx_perm_length(g) = sum(thisgene);
  end
  ridx = 1;
  for g1=1:ngc-1, for g2=g1+1:ngc
    m1 = idx_perm{g1};
    m2 = idx_perm{g2};
    pc = intersect(Tp(m1,2),Tp(m2,2));
    Xc = sum(ismember(Tp(m1,2),pc)) + sum(ismember(Tp(m2,2),pc));
    Xe = idx_perm_length(g1)+idx_perm_length(g2)-Xc;
    R.c(ridx) = R.c(ridx) + (Xc>=R.Nc(ridx) & Xe>=R.Ne(ridx));
    R.pc_ct(ridx) = R.pc_ct(ridx) + (length(pc) >= R.len_pc_obs(ridx));
    ridx=ridx+1;
  end,end
end

% summarize results

R.p = R.c / nruns;
R.pc_frac = R.pc_ct / nruns;
R.relationship = repmat({'none'},npairs,1);
idx = find(R.p<=P.correlation_significance_threshold);
R.relationship(idx) = repmat({'anti-correlated'},length(idx),1);
idx2 = idx(R.pc_frac(idx)<0.5);
R.relationship(idx2) = repmat({'correlated'},length(idx2),1);
R.q = calc_fdr_value(R.p);
R = sort_struct(R,{'q','p'});

if ~isempty(P.gene_gene_correlations_report_filename)
  save_struct(R,P.gene_gene_correlations_report_filename);
end


return

%%%

i=1;
fprintf('%s vs %s\n', R.gene1{i},R.gene2{i});
idx1 = find(strcmp(M.mut.gene_name,R.gene1{i}));
idx2 = find(strcmp(M.mut.gene_name,R.gene2{i}));
sortrows([M.mut.patient_name(idx1) M.mut.Start_position(idx1) M.mut.type(idx1) M.mut.i_tumor_f(idx1)])
sortrows([M.mut.patient_name(idx2) M.mut.Start_position(idx2) M.mut.type(idx2) M.mut.i_tumor_f(idx2)])
