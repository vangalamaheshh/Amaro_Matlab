function M = analyze_silent_nonsilent_ratios(M,P)
% M = analyze_silent_nonsilent_ratios(M,P)
%
% based on analyze_KaKs_data2.m (method 2)
%
% for each gene,
%   given:
%      (1) the observed total number of mutations in each point-mut category (from M.n_nonsilent_ignoring_null_categ, M.n_silent)
%      (2) the expected ratio of nonsilent/silent in each point-mut category (from M.N_non_cov ./ M.N_sil_cov)
%
%   does:
%      (1) in each category, distributes observed mutation counts between the effect categories
%      (2) convolutes distributions to build distribution of NS count
%      (3) trivially derives distribution of S count, and NS/S ratio
%   computes:
%      (1) the expected (most-likely) NS/S ratio
%      (2) the observed NS/S ratio (may include indels or not, depending on setting of parameter)
%      (3) the right-sided p-value (probability of getting at least this high of a NS/S ratio by chance)
%
% Mike Lawrence 2011


% examples:
%
% (1) 40 nonsilent, 0 silent
%     --> strong evidence of unexpectedly high NS/S ratio
% (2) 2 nonsilent, 0 silent
%     --> weak evidence of elevated NS/S ratio
% (3) 40 nonsilent, 15 silent
%     --> strong evidence of expected NS/S ratio
% (4) 3 nonsilent, 1 silent
%     --> weak evidence of expected NS/S ratio

% ask two questions:
%   Q1:  What is the probability that the NS/S ratio is higher than due to chance?
%   Q2:  What is the probability that the NS/S ratio is exactly what we expect by chance?
%   (haven't fully explored these)


if ~exist('P','var'), P=[]; end

P=impose_default_value(P,'use_sample_specific_mutation_rates',false);
P=impose_default_value(P,'impute_full_coverage',~isfield(M,'N_sil_cov'));
P=impose_default_value(P,'pval_cutoff',1e-11);
P=impose_default_value(P,'method',1);
P=impose_default_value(P,'include_indels_in_observed_ns_s_ratio',P.method==2);

if isfield(P,'genes_to_analyze') && ~isempty(P.genes_to_analyze)
  gta = P.genes_to_analyze;
  if ischar(gta), gta={gta}; end
  if islogical(gta), gta=find(gta); end
  if iscellstr(gta)
      gta = listmap(gta,M.gene.name);
  elseif isnumeric(gta)
    % OK
  else
    error('unknown format for P.genes_to_analyze');
  end
else
  gta = 1:M.ng;
end

%if P.use_sample_specific_mutation_rates
%  fprintf('Warning: "use_sample_specific_mutation_rates" not yet implemented for NS/S ratio\n');
%end

fprintf('Analyzing silent/nonsilent ratios\n');

demand_fields(M,{'n_nonsilent_ignoring_null_categ','n_silent'});
if ~P.impute_full_coverage
  demand_fields(M,{'N_sil_cov','N_non_cov'});
else
  demand_fields(M,{'N_sil_terr','N_non_terr'});
end

ncat_including_indels = M.TOT - 1;
ncat = ncat_including_indels - M.NUM_INDEL_CLASSES;
if P.include_indels_in_observed_ns_s_ratio
  ncat_to_include_in_ratio_obs = ncat_including_indels;
else
  ncat_to_include_in_ratio_obs = ncat;
end  

% default results
M.gene.nnon_obs = zeros(M.ng,1);
M.gene.nsil_obs = zeros(M.ng,1);
M.gene.ratio_ns_s_exp = nan(M.ng,1);
M.gene.ratio_ns_s_obs = nan(M.ng,1);
M.gene.pval_ns_s = ones(M.ng,1);

switch P.method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2

for g=gta, if ~mod(g,1000), fprintf('%d/%d ',g,M.ng); end
  nnon = sum(M.n_nonsilent_ignoring_null_categ(g,1:ncat_to_include_in_ratio_obs,:),3);
  nsil = sum(M.n_silent(g,1:ncat_to_include_in_ratio_obs,:),3);

  % calculate expected non/sil ratio
  %     by taking weighted mean of category-specific non/sil expected ratios (weights = observed mutations)

  if ~P.impute_full_coverage
    Nnon = sum(M.N_non_cov(g,1:ncat_to_include_in_ratio_obs,:),3);
    Nsil = sum(M.N_sil_cov(g,1:ncat_to_include_in_ratio_obs,:),3);
  else
    Nnon = M.N_non_terr(g,1:ncat_to_include_in_ratio_obs).*M.np;
    Nsil = M.N_sil_terr(g,1:ncat_to_include_in_ratio_obs).*M.np;
  end
  fnon_exp = Nnon./(Nnon+Nsil);
  indelnull = grepi('indel|null',M.categ.name(1:ncat_to_include_in_ratio_obs),1);
  fnon_exp(indelnull) = 1;
  fnon_exp = weighted_mean(fnon_exp,nnon);
  
  % calculated observed ratio and confidence interval
  nnon = sum(nnon);
  nsil = sum(nsil);
  [fnon_obs fnon_obs_ci] = binofit(nnon,nnon+nsil);
  fnon_obs_stdev = max((fnon_obs-fnon_obs_ci(1)),(fnon_obs_ci(2)-fnon_obs))/1.96;

  % calculate Z-score of observed ratio
  z = (fnon_obs - fnon_exp) / fnon_obs_stdev;  

  if 0
  pval = normcdf(fnon_exp,fnon_obs,fnon_obs_stdev);


  Z = [];  Pval=[];
  for nnon=0:200
    for nsil=0:200
      [fnon_obs fnon_obs_ci] = binofit(nnon,nnon+nsil);
      fnon_obs_stdev = max((fnon_obs-fnon_obs_ci(1)),(fnon_obs_ci(2)-fnon_obs))/1.96;
      z = (fnon_obs - fnon_exp) / fnon_obs_stdev;
      pval = normcdf(fnon_exp,fnon_obs,fnon_obs_stdev);
      Z(nnon+1,nsil+1) = z;
      Pval(nnon+1,nsil+1) = pval;
  end,end
  end



end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1

for g=gta, if ~mod(g,1000), fprintf('%d/%d ',g,M.ng); end
  nnon = sum(M.n_nonsilent_ignoring_null_categ(g,:,:),3);
  nsil = sum(M.n_silent(g,:,:),3);
  ntot = nnon+nsil;
  if sum(ntot(1:ncat))>0   % gene has at least one point mutation (not counting indels)
    if ~P.impute_full_coverage
      Nnon = sum(M.N_non_cov(g,:,:),3);
      Nsil = sum(M.N_sil_cov(g,:,:),3);
    else
      Nnon = M.N_non_terr(g,:).*M.np;
      Nsil = M.N_sil_terr(g,:).*M.np;
    end
    Ntot = Nnon + Nsil;
    if any(Ntot==0 & ntot>0)
      fprintf('WARNING: failed to remove impossible mutation in %s\n',M.gene.name{g});
      continue;
    end
    fNS = Nnon ./ Ntot;
    NS = {};
    for c=1:ncat
      if ntot(c)>0
        if fNS(c)>=0 && fNS(c)<=1
          NS{end+1,1} = binopdf(0:ntot(c),ntot(c),fNS(c));
        else
          fprintf('Unexpected fNS = %f in analyze_silent_nonsilent_ratios\n',fNS);
          keyboard
        end
      end
    end
    % convolute to get final distribution
    dist = [];
    dist.prob = batch_convolute(NS);
    dist.nnon = (0:slength(dist)-1)';
    dist.nsil = slength(dist)-1-dist.nnon;
    dist.ratio = dist.nnon ./ dist.nsil;

    % MLE
    [tmp idx] = max(dist.prob);
    nnon_exp = dist.nnon(idx);
    nsil_exp = dist.nsil(idx);
    ratio_exp = dist.ratio(idx);

    % p-value
    nnon_obs_tot = sum(nnon(1:ncat_to_include_in_ratio_obs));
    nsil_obs_tot = sum(nsil(1:ncat_to_include_in_ratio_obs));
    ratio_obs = nnon_obs_tot / nsil_obs_tot;   
    idx = find(dist.ratio>=ratio_obs);
    if isempty(idx)
      fprintf('Unexpected case in analyze_silent_nonsilent_ratios\n');
      keyboard
      pval = 0;
    elseif idx==1
      pval = 1;
    else
      pval = 1-sum(dist.prob(1:idx-1));     
    end

    % record results
    M.gene.nnon_obs(g) = nnon_obs_tot;
    M.gene.nsil_obs(g) = nsil_obs_tot;
    M.gene.ratio_ns_s_obs(g) = ratio_obs;
    M.gene.ratio_ns_s_exp(g) = ratio_exp;
    M.gene.pval_ns_s(g) = pval;
    
    % pause to show plot for genes-of-interest
    goi = {};
    %goi = {'BCL2'};
    %  goi = {'TP53','PPFIA1','NFE2L2','TTN','CDH10','PLA2G2A','TRIM39','SI','FLG',...
    %       'RYR2','CSMD3','ZNF479','NLRP8','TOLLIP','EPHB3','EPHB1','CDKN2A'};
    %  goi = {'KEAP1','MUC17','DNAH5','NOTCH1','PIK3CA','PIK3CG','COL11A1'};
    if ismember(M.gene.name{g},goi)
      figure(1);clf;for i=1:length(NS)
        subplot(length(NS),1,i);bar(0:length(NS{i})-1,NS{i});xlim([-0.5 length(NS{i})-0.5]);
        text(nnon(i),NS{i}(nnon(i)+1)+0.1,['obs=' num2str(nnon(i))],'color',[1 0 0],'horizontalalignment','center');
        title(M.mutclass{i});end;set(gcf,'color',[1 1 1]);
      figure(2)
      bar(dist.nnon,dist.prob),xlim([-0.5 round(max(dist.nnon)*1.3)]),colormap([0.7 0.7 0.7]);
      pmax = max(dist.prob);
      title(M.gene.name{g},'fontsize',30);
      xlabel('# NS','fontsize',20); ylabel('prob','fontsize',20);
      line((nnon_obs_tot-0.5) * [1 1],[0 pmax*1.2],'color',[1 0 0],'linewidth',2);
      text(nnon_obs_tot * 1.02,pmax*1.06,sprintf('p = %0.2f',pval),'color',[1 0 0],'fontsize',20);
      text(nnon_obs_tot * 0.2,pmax*1.06,sprintf('nonsilent = %d',nnon_obs_tot),'fontsize',12);
      text(nnon_obs_tot * 0.2,pmax*1.00,sprintf('silent = %d',nsil_obs_tot),'fontsize',12);
      text(nnon_obs_tot * 0.2,pmax*0.94,sprintf('total = %d',nnon_obs_tot+nsil_obs_tot),'fontsize',12);
      xlim([-0.5 nnon_obs_tot+nsil_obs_tot+0.5]);
      ylim([0 pmax*1.2]);
      keyboard
    end
  end  % if ntot>0

end  % next gene

end
%%%%%%%%%%%% switch P.method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% fix negative and unreliable very-small values
M.gene.pval_ns_s_lessthan_flag = (M.gene.pval_ns_s <= P.pval_cutoff);
M.gene.pval_ns_s(M.gene.pval_ns_s_lessthan_flag) = P.pval_cutoff;

fprintf('Done\n');

