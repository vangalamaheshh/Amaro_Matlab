function R = analyze_gene_gene_correlations(M,P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'genes_to_consider','*required*');
P = impose_default_value(P,'restrict_to_one_mut_per_gene',true);
P = impose_default_value(P,'gene_gene_correlations_num_permutations',10000);
P = impose_default_value(P,'gene_gene_correlations_report_filename',[]);
P = impose_default_value(P,'correlation_significance_threshold',0.05);
P = impose_default_value(P,'max_ci_ratio',0.5);
P = impose_default_value(P,'randseed',1234);

demand_fields(M,{'mut','gene','use_nonsilent'});
demand_fields(M.mut,{'gene','patient'});
demand_fields(M.gene,'name');

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

% compute observed statistics

npairs = (ngc * (ngc-1)) / 2;

T = [M.mut.gene(M.use_nonsilent) M.mut.patient(M.use_nonsilent)];
if P.restrict_to_one_mut_per_gene, T = unique(T,'rows'); end

maxpt = max(T(:,2));
patct = histc(T(:,2),1:maxpt);

maxgene = max(T(:,1));
genect = histc(T(:,1),1:maxgene);

idx_obs = cell(ngc,1); for g=1:ngc, idx_obs{g} = find(T(:,1)==gc(g)); end
R=[]; ridx = 0;
for g1=1:ngc-1, for g2=g1+1:ngc,ridx=ridx+1;
  R.gene1idx(ridx,1) = g1;
  R.gene2idx(ridx,1) = g2;
  R.gene1{ridx,1} = M.gene.name{gc(g1)};
  R.gene2{ridx,1} = M.gene.name{gc(g2)};
  R.gene1ct(ridx,1) = genect(gc(g1));
  R.gene2ct(ridx,1) = genect(gc(g2));
  m1_obs = idx_obs{g1};
  m2_obs = idx_obs{g2};
  pats_with_both = intersect(T(m1_obs,2),T(m2_obs,2));
  R.npats_with_both(ridx,1) = length(pats_with_both);
  R.nmuts_together(ridx,1) = sum(ismember(T([m1_obs;m2_obs],2),pats_with_both));
  R.nmuts_separate(ridx,1) = length([m1_obs;m2_obs])-R.nmuts_together(ridx);
  R.nperms(ridx,1) = 0;
  R.nperms_asgood(ridx,1) = 0;
  R.nperms_together(ridx,1) = 0;
  R.done(ridx,1) = false;
end,end

% do permutations

fprintf('Performing permutations:\n');
rand('twister',P.randseed);
nperms = P.gene_gene_correlations_num_permutations;

pats = cell(ngc,1); for g=1:ngc, pats{g} = nan(genect(gc(g)),1); end

ndone = 0;
ndone_gene = zeros(ngc,1);
reportrate = 1; % seconds
tt=tic; oldto = 0; oldreport = 0; stepsize = 10; nextreport = stepsize;
for i=1:nperms

  if i==nextreport || ndone==ngc
    to = toc(tt);
    stepwas = (i-oldreport);
    fprintf('%.2f\t%d/%d\t%.2f\t%.3f\t%.1f\t%d/%d pairs\t%d/%d genes done\n',to,i,nperms,to/i,...
       (to-oldto)/stepwas,stepwas/(to-oldto),...
       ndone,npairs,sum(ndone_gene==ngc-1),ngc);
    if (to-oldto)<reportrate, stepsize = round(stepsize * 1.2); end
    oldto = to;
    oldreport = i;
    nextreport = nextreport + stepsize;
  end
  
  % generate permutation

  gidx = find(ndone_gene<ngc-1);
  if isempty(gidx), break; end
  for gi=1:length(gidx), g=gidx(gi);
    xx = [(1:maxpt)' patct];
    for p=1:genect(gc(g))
      xx(:,3) = cumsum(xx(:,2));
      r = ceil(rand*xx(end,3));
      row = find(xx(:,3)>=r,1,'first');
      pats{g}(p) = xx(row,1);
      xx(row,2) = 0;
    end
  end

  % analyze pairs

  ridx = find(~R.done);
  for ri=1:length(ridx), r=ridx(ri);

    g1 = R.gene1idx(r);
    g2 = R.gene2idx(r);
    pat1 = pats{g1};
    pat2 = pats{g2};

    [pats_with_both pci pcj] = intersect(pat1,pat2);
    npats_with_both = length(pats_with_both);
    nmuts_together = length(pci)*2;
    nmuts_separate = genect(gc(g1)) + genect(gc(g2)) - nmuts_together;
    is_asgood = (nmuts_together>=R.nmuts_together(r) && nmuts_separate>=R.nmuts_separate(r));
%    fprintf('nt_obs %d    ns_obs %d   nt %d    ns %d    is_asgood %d\n',...
%            R.nmuts_together(r),R.nmuts_separate(r),nmuts_together,nmuts_separate,is_asgood);
    if is_asgood && (nmuts_together>R.nmuts_together(r) ||  nmuts_separate>R.nmuts_separate(r))
      keyboard
    end
    R.nperms_asgood(r) = R.nperms_asgood(r) + is_asgood;
    R.nperms_together(r) = R.nperms_together(r) + (npats_with_both >= R.npats_with_both(r));

    ci_ratio = calc_pval_ci_ratio(R.nperms_asgood(r),i);
    if ci_ratio <= P.max_ci_ratio
      R.nperms(r) = i;
      R.done(r) = true;
      ndone = ndone + 1;
      ndone_gene(g1) = ndone_gene(g1) + 1;
      ndone_gene(g2) = ndone_gene(g2) + 1;
    end
  end

end

R.nperms(~R.done) = nperms;

% summarize results

R.p = R.nperms_asgood ./ R.nperms;
R.frac_perms_together = R.nperms_together ./ R.nperms;
R.relationship = repmat({'none'},npairs,1);
idx = find(R.p<=P.correlation_significance_threshold);
R.relationship(idx) = repmat({'anti-correlated'},length(idx),1);
idx2 = idx(R.frac_perms_together(idx)<0.5);
R.relationship(idx2) = repmat({'correlated'},length(idx2),1);
R.q = calc_fdr_value(R.p);
R = sort_struct(R,{'q','p'});

if ~isempty(P.gene_gene_correlations_report_filename)
  save_struct(R,P.gene_gene_correlations_report_filename);
end

idx=1:min(npairs,30);
[R.gene1(idx) R.gene2(idx) num2cell([R.p(idx) R.q(idx)]) R.relationship(idx)]


return

%%%

i=1;
fprintf('%s vs %s\n', R.gene1{i},R.gene2{i});
idx1 = find(strcmp(M.mut.gene_name,R.gene1{i}));
idx2 = find(strcmp(M.mut.gene_name,R.gene2{i}));
sortrows([M.mut.patient_name(idx1) M.mut.Start_position(idx1) M.mut.type(idx1) M.mut.i_tumor_f(idx1)])
sortrows([M.mut.patient_name(idx2) M.mut.Start_position(idx2) M.mut.type(idx2) M.mut.i_tumor_f(idx2)])
