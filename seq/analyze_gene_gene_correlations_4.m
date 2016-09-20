function R = analyze_gene_gene_correlations_4(M,P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'genes_to_consider','*required*');
P = impose_default_value(P,'gene_gene_correlations_num_permutations',1e6);
P = impose_default_value(P,'gene_gene_correlations_report_filename',[]);
P = impose_default_value(P,'gene_gene_correlations_verbose_report',true);
P = impose_default_value(P,'correlation_significance_threshold',0.05);
P = impose_default_value(P,'max_ci_ratio',2);
P = impose_default_value(P,'randseed',1234);
P = impose_default_value(P,'method',2);
P = impose_default_value(P,'debug',false);

demand_fields(M,{'mut','gene'});
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
  R=[];
  fprintf('Need at least 2 genes to consider.\n');
  return
end

fprintf('Analyzing gene-gene correlations and anti-correlations\n');

% compute observed statistics

npairs = (ngc * (ngc-1)) / 2;

if isfield(M,'use_nonsilent')
  useidx = M.use_nonsilent;
elseif isfield(M.mut,'is_coding') && isfield(M.mut,'is_silent')
  useidx = find(M.mut.is_coding & ~M.mut.is_silent);
elseif isfield(M.mut,'type')
  useidx = find(~ismember(M.mut.type,{'Synonymous','Silent','Intron','UTR','3''-UTR','5''-UTR'}));
else
  fprintf('Unable to find information to distinguish nonsilent coding mutations... using all mutations.\n');
  useidx = (1:slength(M.mut))';
end

if isempty(useidx), error('No mutations!'); end

if isnumeric(M.mut.gene) && isnumeric(M.mut.patient)
  T = unique([M.mut.gene(useidx) M.mut.patient(useidx)],'rows');
else
  T = unique([M.mut.gene_idx(useidx) M.mut.pat_idx(useidx)],'rows');
end

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

% Fisher's exact test
nboth = R.npats_with_both;
nonly1 = R.gene1ct - nboth;
nonly2 = R.gene2ct - nboth;
neither = maxpt - nboth - nonly1 - nonly2;
R.pval_fisher = nan(npairs,1);
for i=1:npairs, R.pval_fisher(i) = fisher_exact_test(nboth(i),nonly1(i),nonly2(i),neither(i)); end
R.sig_fisher = (R.pval_fisher<=P.correlation_significance_threshold);
R.relationship_fisher = repmat({'none'},npairs,1);
idx = find(R.sig_fisher);
R.relationship_fisher(idx) = repmat({'anti-correlated'},length(idx),1);
idxa = intersect(idx,find(nboth>nonly1.*nonly2./neither));
R.relationship_fisher(idxa) = repmat({'correlated'},length(idxa),1);

% permutations
R.nperms = zeros(npairs,1);
R.nperms_as_corr = zeros(npairs,1);
R.nperms_as_anti = zeros(npairs,1);
R.corr_done = false(npairs,1);
R.anti_done = false(npairs,1);
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

fprintf('sec  perm/nperms  sec/perm   sec/perm(local)  perm/sec(local)\n');

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
    R.nperms_as_corr(r) = R.nperms_as_corr(r) + (nmuts_together>=R.nmuts_together(r));
    R.nperms_as_anti(r) = R.nperms_as_anti(r) + (nmuts_separate>=R.nmuts_separate(r));

    R.corr_done(r) = (calc_pval_ci_ratio(R.nperms_as_corr(r),i) <= P.max_ci_ratio);
    R.anti_done(r) = (calc_pval_ci_ratio(R.nperms_as_anti(r),i) <= P.max_ci_ratio);
    if ~R.sig_fisher(r) || (R.corr_done(r) && R.anti_done(r))
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
R.pval_perms_corr = (1+R.nperms_as_corr) ./ (2+R.nperms);
R.pval_perms_corr_lessthan_flag = (R.nperms_as_corr==0);
R.pval_perms_anti = (1+R.nperms_as_anti) ./ (2+R.nperms);
R.pval_perms_anti_lessthan_flag = (R.nperms_as_anti==0);

R.relationship = repmat({'none'},npairs,1);
R.pval = R.pval_fisher;
R.pval_lessthan_flag = false(npairs,1);

idx = find(strcmp('correlated',R.relationship_fisher));
R.pval(idx) = max(R.pval_fisher(idx),R.pval_perms_corr(idx));
idx2 = idx(R.pval(idx)==R.pval_perms_corr(idx));
R.pval_lessthan_flag(idx2) = R.pval_perms_corr_lessthan_flag(idx2);
idx2 = idx(R.pval(idx)<=P.correlation_significance_threshold);
R.relationship(idx2) = repmat({'correlated'},length(idx2),1);

idx = find(strcmp('anti-correlated',R.relationship_fisher));
R.pval(idx) = max(R.pval_fisher(idx),R.pval_perms_anti(idx));
idx2 = idx(R.pval(idx)==R.pval_perms_anti(idx));
R.pval_lessthan_flag(idx2) = R.pval_perms_anti_lessthan_flag(idx2);
idx2 = idx(R.pval(idx)<=P.correlation_significance_threshold);
R.relationship(idx2) = repmat({'anti-correlated'},length(idx2),1);

R.qval = calc_fdr_value(R.pval);
R = sort_struct(R,{'qval','pval'});

fprintf('\nTop gene-gene pairs:\n');
fprintf('gene1\tgene2\tp\tq\trelationship\n');
maxtoshow = 60;
for i=1:min(npairs,maxtoshow)
  if R.pval_lessthan_flag(i), lt='<'; else lt=''; end
  fprintf('%s\t%s\t%s%f\t%s%f\t%s\n',...
          R.gene1{i}, R.gene2{i}, lt, R.pval(i), lt, R.qval(i), R.relationship{i});
end

if ~P.gene_gene_correlations_verbose_report
  R = rename_fields(R,{'pval','qval','pval_lessthan_flag'},{'p','q','p_lessthan_flag'});
  R = keep_fields(R,{'gene1','gene2','p_lessthan_flag','p','q','relationship'});
end

if ~isempty(P.gene_gene_correlations_report_filename)
  save_struct(R,P.gene_gene_correlations_report_filename);
end
