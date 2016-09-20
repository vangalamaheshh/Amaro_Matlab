function R = analyze_gene_gene_correlations(M,P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'genes_to_consider','*required*');
P = impose_default_value(P,'gene_gene_correlations_num_permutations',10000);
P = impose_default_value(P,'gene_gene_correlations_report_filename',[]);
P = impose_default_value(P,'correlation_significance_threshold',0.05);
P = impose_default_value(P,'max_ci_ratio',0.5);
P = impose_default_value(P,'randseed',1234);
P = impose_default_value(P,'method',2);
P = impose_default_value(P,'debug',false);

demand_fields(M,{'mut','gene','use_nonsilent'});
demand_fields(M.mut,{'gene','patient'});
demand_fields(M.gene,'name');

gc = listmap(P.genes_to_consider,M.gene.name);
notfound = isnan(gc);
if any(notfound)
  fprintf('Genes not found:');disp(gc(notfound));
end
gc = gc(~notfound);
gcname = P.genes_to_consider(~notfound);
ngc = length(gc);
if ngc<2
  fprintf('Need at least 2 genes to consider.\n');
  return
end

% compute observed statistics

npairs = (ngc * (ngc-1)) / 2;

T = unique([M.mut.gene(M.use_nonsilent) M.mut.patient(M.use_nonsilent)],'rows');
%   (restrict to one mutation per gene per patient)

% M no longer needed in rest of function

maxpt = max(T(:,2));
patct = histc(T(:,2),1:maxpt);

maxgene = max(T(:,1));
genect = histc(T(:,1),1:maxgene);
idx_obs = cell(ngc,1); for g=1:ngc, idx_obs{g} = find(T(:,1)==gc(g)); end
obs_pats = cell(ngc,1); for g=1:ngc, obs_pats{g} = T(idx_obs{g},2); end

R=[]; ridx = 0;
R.gene1idx = nan(npairs,1);
R.gene2idx = nan(npairs,1);
R.gene1 = [];
R.gene2 = [];
R.gene1ct = [];
R.gene2ct = [];
R.npats_with_both = nan(npairs,1);
for g1=1:ngc-1, for g2=g1+1:ngc,ridx=ridx+1;
  R.gene1idx(ridx) = g1;
  R.gene2idx(ridx) = g2;
  m1_obs = idx_obs{g1};
  m2_obs = idx_obs{g2};
  R.npats_with_both(ridx) = length(intersect(T(m1_obs,2),T(m2_obs,2)));
end,end
R.gene1 = gcname(R.gene1idx);
R.gene2 = gcname(R.gene2idx);
R.gene1ct = genect(gc(R.gene1idx));
R.gene2ct = genect(gc(R.gene2idx));
R.nmuts_together = R.npats_with_both*2;
R.nmuts_separate = R.gene1ct+R.gene2ct-R.nmuts_together;
R.nperms = zeros(npairs,1);
R.nperms_as_together = zeros(npairs,1);
R.nperms_as_separate = zeros(npairs,1);
R.together_done = false(npairs,1);
R.separate_done = false(npairs,1);
R.done = false(npairs,1);

if P.debug
  xobs = zeros(ngc,maxpt); for g=1:ngc, xobs(g,obs_pats{g}) = 1; end
  [tmp xord] = sort(patct);
  clf
  subplot(6,1,1);
    bar(patct(xord)); xlim([1 maxpt]);
    set(gca,'tickdir','out','xtick',[],'ytick',[]);
    ylabel('mutrate');
  subplot(6,1,2);
    imagesc(xobs(:,xord)); xlim([1 maxpt]);
    set(gca,'tickdir','out','ytick',1:ngc,'yticklabel',gcname,'xtick',[]);
    title('observed','fontsize',10);
end

% do permutations

fprintf('Performing permutations:\n');
rand('twister',P.randseed);
nperms = P.gene_gene_correlations_num_permutations;

pats = cell(ngc,1); for g=1:ngc, pats{g} = nan(genect(gc(g)),1); end

ndone = 0;
ndone_gene = zeros(ngc,1);
reportrate = 5; % seconds
tt=tic; oldto = 0; oldreport = 0; stepsize = 10; nextreport = stepsize;
for i=1:nperms

  if i>=nextreport || ndone==npairs
    to = toc(tt);
    stepwas = (i-oldreport);
    fprintf('%.2f\t%d/%d\t%.2f\t%.3f\t%.1f\t%d/%d pairs\t%d/%d genes done\n',to,i,nperms,to/i,...
       (to-oldto)/stepwas,stepwas/(to-oldto),...
       ndone,npairs,sum(ndone_gene==ngc-1),ngc);
    if (to-oldto)<reportrate, stepsize = stepsize * reportrate / (to-oldto); end
    oldto = to;
    oldreport = i;
    nextreport = nextreport + stepsize;
  end
  
  if ndone==npairs, break; end  % finished

  % generate permutation

  if P.method==2
  %     for each gene, randomize the patient labels by weighted random permutation without replacement
    gidx = find(ndone_gene<ngc-1);
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

  elseif P.method==1   % original TSP paper method
  %     randomly permute patient labels (allows two mutations to end up in same gene+patient)    
    Tp = [T(randperm(size(T,1)),1) T(:,2)];
    for g=1:ngc, pats{g} = Tp(Tp(:,1)==gc(g),2); end

  elseif P.method==3   % Fisher's exact test
    nboth = R.npats_with_both;
    nonly1 = R.gene1ct - nboth;
    nonly2 = R.gene2ct - nboth;
    neither = maxpt - nboth - nonly1 - nonly2;
    fishp = nan(npairs,1);
    for i=1:npairs, fishp(i) = fisher_exact_test(nboth(i),nonly1(i),nonly2(i),neither(i)); end
    R.nperms = ones(npairs,1);
    R.nperms_as_together = ones(npairs,1);
    R.nperms_as_separate = ones(npairs,1);
    idx = find((nboth+neither)>(nonly1+nonly2));
    R.nperms_as_together(idx) = fishp(idx);
    idx = find((nboth+neither)<=(nonly1+nonly2));
    R.nperms_as_separate(idx) = fishp(idx);
    R.done = true(npairs,1);
    break;
  else
    error('unknown P.method = %s',P.method);
  end

  % debug mode: plot the first 4 permutations
  if P.debug && i<=4
    xperm = zeros(ngc,maxpt); for g=1:ngc, xperm(g,pats{g}) = 1; end
    subplot(6,1,2+i);
      imagesc(xperm(:,xord)); xlim([1 maxpt]);
      set(gca,'tickdir','out','ytick',1:ngc,'yticklabel',gcname,'xtick',[]);
      title(['permutation #' num2str(i)],'fontsize',10);
    if i==4, keyboard; end
  end

  % analyze pairs

  ridx = find(~R.done);
  for ri=1:length(ridx), r=ridx(ri);

    g1 = R.gene1idx(r);
    g2 = R.gene2idx(r);
    pat1 = pats{g1};
    pat2 = pats{g2};

    npats_with_both = length(intersect(pat1,pat2));
    nmuts_together = npats_with_both*2;
    nmuts_separate = genect(gc(g1)) + genect(gc(g2)) - nmuts_together;
    R.nperms_as_together(r) = R.nperms_as_together(r) + (nmuts_together>=R.nmuts_together(r));
    R.nperms_as_separate(r) = R.nperms_as_separate(r) + (nmuts_separate>=R.nmuts_separate(r));

    R.together_done(r) = (calc_pval_ci_ratio(R.nperms_as_together(r),i) <= P.max_ci_ratio);
    R.separate_done(r) = (calc_pval_ci_ratio(R.nperms_as_separate(r),i) <= P.max_ci_ratio);
    if R.together_done(r) && R.separate_done(r)
      R.done(r) = true;
      R.nperms(r) = i;
      ndone = ndone + 1;
      ndone_gene(g1) = ndone_gene(g1) + 1;
      ndone_gene(g2) = ndone_gene(g2) + 1;
    end
  end

end

R.nperms(~R.done) = nperms;

% summarize results

R.p_together = R.nperms_as_together ./ R.nperms;
R.p_separate = R.nperms_as_separate ./ R.nperms;
R.p = min(R.p_together,R.p_separate);
R.relationship = repmat({'none'},npairs,1);
idx = find(R.p<=P.correlation_significance_threshold);
R.relationship(idx) = repmat({'anti-correlated'},length(idx),1);
idx2 = idx(R.p_together(idx)<R.p_separate(idx));
R.relationship(idx2) = repmat({'correlated'},length(idx2),1);
R.q = calc_fdr_value(R.p);
R = sort_struct(R,{'q','p'});

if ~isempty(P.gene_gene_correlations_report_filename)
  save_struct(R,P.gene_gene_correlations_report_filename);
end

idx=1:min(npairs,30);
[R.gene1(idx) R.gene2(idx) num2cell([R.p(idx) R.q(idx)]) R.relationship(idx)]

