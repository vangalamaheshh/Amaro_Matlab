function [G m v] = analyze_mutrate_covariates(G,V,vidx,P)
% [G m v] = analyze_mutrate_covariates(G,V,vidx,P)
%
% analyzes covariates to find Ffit, the relative mutation rate of each gene
%
% given:
%   G = genes (nfit and Nfit required)
%   V = covariates (val required)
%   vidx = which covariates to use
%
% returns
%   G = genes (with Ffit and other fields added)
%   m = median of Ffit
%   v = variance of Ffit
%
% based on /xchip/cga1/lawrence/mut/analysis/20110909_pancan/run.m

if nargin<1, error('need at least G'); end

if nargin==2 && isstruct(V) && ~isfield(V,'val')
  P=V;
  V=[];
end

if nargin==3 && isstruct(vidx)
  P=vidx;
  clear vidx;
end

if ~exist('V','var') || isempty(V)
  V.val = cell(0,1);
  V.divs = cell(0,1);
end

if ~isstruct(G), error('G should be a struct'); end
if ~isstruct(V), error('V should be a struct'); end

% allow case of only one covariate and no cell structure
flds = fieldnames(V);
flds1 = getfield(V,flds{1});
if ~iscell(flds1)
  tmp = V;
  V = [];
  for i=1:length(flds)
    V = setfield(V,flds{i},{getfield(tmp,flds{i})});
  end
end

if ~exist('vidx','var'), vidx = 1:slength(V); end
if any(vidx<1) || any(vidx>slength(V)), error('vidx out of range'); end
V = reorder_struct(V,vidx);


ng = slength(G);
nv = slength(V);

if ~exist('P','var'), P=[];end
P=impose_default_value(P,'method','nearest_neighbors');
P=impose_default_value(P,'generate_plot',false);
P=impose_default_value(P,'genes_to_analyze',1:ng);

if iscellstr(P.genes_to_analyze)
  demand_field(G,'name');
  P.genes_to_analyze = listmap(P.genes_to_analyze,G.name);
end
analyze_this_gene = ismember(1:ng,P.genes_to_analyze);

if nv==0, P.method = 'binning'; end

if strcmpi(P.method,'binning')

  demand_fields(G,{'nfit','Nfit'});
  demand_fields(V,{'val','divs'});

  % first, do binning for each covariate
  V.binmin = cell(nv,1); V.binmax = cell(nv,1); V.nbins = nan(nv,1);
  for vi=1:nv
    V.binmin{vi} = [-inf V.divs{vi}]';
    V.binmax{vi} = [V.divs{vi} +inf]';
    V.nbins(vi) = length(V.binmin{vi});
    for bi=1:V.nbins(vi)
      bin_gidx{vi}{bi} = find(V.val{vi}>=V.binmin{vi}(bi) & V.val{vi}<V.binmax{vi}(bi));
    end
  end

  tot_nbins = prod(V.nbins);
  
  bin=[];
  bin.bini = nan(tot_nbins,nv);
  bin.ct = nan(tot_nbins,1);
  bin.nfit = nan(tot_nbins,1);
  bin.Nfit = nan(tot_nbins,1);
  G.bin = nan(ng,1);

  % now bin the genes
  bini = ones(1,nv);
  for tbi=1:tot_nbins
    % find genes that belong to this bin
    gidx = 1:ng;
    for vi=1:nv
      gidx = intersect(gidx, bin_gidx{vi}{bini(vi)});
    end
    G.bin(gidx) = tbi;
    bin.bini(tbi,:) = bini;
    bin.ct(tbi,1) = length(gidx);
    bin.nfit(tbi,1) = nansum(G.nfit(gidx));
    bin.Nfit(tbi,1) = round(nansum(G.Nfit(gidx)));
    % increment bin
    vi=1;
    if tbi<tot_nbins
      while(true)
        bini(vi)=bini(vi)+1;
        if bini(vi)<=V.nbins(vi), break; end
        bini(vi)=1;
        vi=vi+1;
        if vi>nv, break; end  % last iteration
      end
    end
  end % next bin
  
  [bin.rate bin.ci_low bin.ci_high] = binofit_3d(bin.nfit,max(bin.nfit*1e4,bin.Nfit));
  bin.rate(bin.Nfit==0)=nan; bin.ci_low(bin.Nfit==0)=nan; bin.ci_high(bin.Nfit==0)=nan;
  
  totrate_fit=nansum(G.nfit)/nansum(G.Nfit);
  bin.relrate = bin.rate/totrate_fit; bin.relci_low = bin.ci_low/totrate_fit; bin.relci_high=bin.ci_high/totrate_fit;
  
  G.Ffit = nansub(bin.relrate,G.bin,1);
  
elseif strcmpi(P.method, 'nearest_neighbors')

  P=impose_default_value(P,'nearest_neighbors_binofit_alpha',0.05);
  P=impose_default_value(P,'nearest_neighbors_max_theta',0.1);
  P=impose_default_value(P,'nearest_neighbors_min_neighbs',10);
  P=impose_default_value(P,'nearest_neighbors_add_stepsize',10);
  
  demand_fields(G,{'nfit','Nfit'});
  demand_fields(V,'val');
  
 % first convert each covariate to Z-scores and replace missing values with Z=0
  Z = nan(ng,nv);
  for vi=1:nv
    V.missing{vi,1} = (isnan(V.val{vi}) | isinf(V.val{vi}));
    mn = mean(V.val{vi}(~V.missing{vi}));
    sd = std(V.val{vi}(~V.missing{vi}));
    Z(:,vi) = (V.val{vi}-mn)./sd;
    Z(V.missing{vi},vi) = 0;
  end

  % now, for each gene, find enough nearest-neighbors to satisfy the requirements
  G.fitrate = nan(ng,1);
  G.fit_ng = nan(ng,1);
  G.fit_n = nan(ng,1);
  G.fit_N = nan(ng,1);
  G.fit_d = nan(ng,1);
  for gi=1:ng, if ~mod(gi,1000), fprintf('%d/%d ',gi,ng); end
    if ~analyze_this_gene(gi), continue; end

    % compute and sort squared-distances
    dist2 = sum(bsxfun(@minus,Z,Z(gi,:)).^2,2);
    [tmp ord] = sort(dist2);
    ord(ord==gi)=[];
    method=2;
    if method==1
      ct = P.nearest_neighbors_min_neighbs;
      n = sum(G.nfit(ord(1:ct)));
      N = sum(G.Nfit(ord(1:ct)));
      while(true)
        theta = calc_pval_ci_ratio(n,N);
        %      [r ci] = binofit(n,N,P.nearest_neighbors_binofit_alpha);
        %      theta = (ci(2)-ci(1))/r;
        if (theta<=P.nearest_neighbors_max_theta), break; end
        if (ct==ng), break; end
        ct=ct+1;
        n=n+G.nfit(ord(ct));
        N=N+G.Nfit(ord(ct));
      end
      G.fitrate(gi) = n/N;
    else % method==2
      n = cumsum(G.nfit(ord));
      N = cumsum(G.Nfit(ord));
      ngenes = (1:length(n))';
      d = dist2(ord);
      keep = P.nearest_neighbors_min_neighbs:P.nearest_neighbors_add_stepsize:length(n);
      n = n(keep);
      N = N(keep);
      ngenes = ngenes(keep);
      d = d(keep);
      tmp = (N-n+1) ./ ((n+1).* (N+3));
      tmp(tmp<0) = inf;
%      theta = 2*1.96*sqrt(tmp);
%      idx = find(theta<=P.nearest_neighbors_max_theta,1);
      idx = find(tmp<=(P.nearest_neighbors_max_theta/(2*1.96)).^2,1);
      if isempty(idx), idx = length(n); end
      G.fit_ng(gi) = ngenes(idx);
      G.fit_n(gi) = n(idx);
      G.fit_N(gi) = N(idx);
      G.fit_d(gi) = d(idx);
      G.fitrate(gi) = n(idx)/N(idx);
    end
    
  end, fprintf('\n'); % next gene

  totrate_fit=nansum(G.nfit)/nansum(G.Nfit);
  G.Ffit = G.fitrate / totrate_fit;

elseif strcmpi(P.method, 'newmethod1_failed_attempt1')

  P=impose_default_value(P,'newmethod1_Frange_to_test','*required*');
  P=impose_default_value(P,'newmethod1_num_neighbors','*required*');
  P=impose_default_value(P,'newmethod1_min_dist2','*required*');

  Frange = P.newmethod1_Frange_to_test;
  mindist2 = P.newmethod1_min_dist2;

  demand_fields(G,{'name','nnon','nsil','Nnon','Nsil'});
  demand_fields(V,'val');

  globalrate_non = sum(G.nnon)/sum(G.Nnon);
  globalrate_sil = sum(G.nsil)/sum(G.Nsil);

 % first convert each covariate to Z-scores and replace missing values with Z=0
  Z = nan(ng,nv);
  for vi=1:nv
    V.missing{vi,1} = (isnan(V.val{vi}) | isinf(V.val{vi}));
    mn = mean(V.val{vi}(~V.missing{vi}));
    sd = std(V.val{vi}(~V.missing{vi}));
    Z(:,vi) = (V.val{vi}-mn)./sd;
    Z(V.missing{vi},vi) = 0;
  end

  % now, for each gene, find distribution on F
  G.Ffit_distrib = cell(ng,1);
  G.details = cell(ng,1);

  for gi=1:ng, if ~mod(gi,1000), fprintf('%d/%d ',gi,ng); end
    if ~analyze_this_gene(gi), continue; end

    % compute and sort squared-distances
    dist2 = sum(bsxfun(@minus,Z,Z(gi,:)).^2,2);
    [tmp ord] = sort(dist2);
    ord(ord==gi)=[]; % remove this gene itself
    dist2 = max(dist2,mindist2);

    % add contribution of each neighbor
    %    -- get likelihood across distribution
    %    -- (could normalize to max, but won't do this for now)
    %    -- multiply by N^2 that the estimate is based on
    %    -- divide by Distance^2 from the gene itself


    Fdist = zeros(size(Frange));
    details = {}; didx = 1;

    for ni=0:P.newmethod1_num_neighbors
      if ni==0
        gidx = gi;
        dist2_g = mindist2;
      else
        gidx = ord(ni);
        dist2_g = dist2(gidx);
      end
      name = G.name{gidx};

      for muttype=1:2

        if muttype==1
          type = 'silent';
          n = G.nsil(gidx);
          N = G.Nsil(gidx);
          globalrate = globalrate_sil;
        elseif muttype==2
          type = 'nonsilent';
          if ni==0, continue; end  % don't count nonsilent mutations of the gene itself
          n = G.nnon(gidx);
          N = G.Nnon(gidx);
          globalrate = globalrate_non;
        end

        if (N>n)
%          if strcmp(name,'KEAP1'), keyboard; end
          lk = binopdf(n,N,globalrate*Frange);
          lk = lk/sum(lk);
          factor = 1; %(N.^2)/dist2_g;
          Fdist = Fdist + lk*factor;
          details.source{didx,1} = type;
          details.ni(didx,1) = ni;
          details.name{didx,1} = name;
          details.dist2(didx,1) = dist2_g;
          details.factor(didx,1) = factor;
          details.n(didx,1) = n;
          details.N(didx,1) = N;
          details.F(didx,1) = (n./N)/globalrate;
          [tmp idx] = max(Fdist);
          details.Fcum(didx,1) = Frange(idx);
          details.Fend(didx,1) = Fdist(end)/sum(Fdist);
          didx=didx+1;
        end
      end % next mutation type
    end % next neighbor
    
    % normalize the distribution
    Fdist = Fdist / sum(Fdist);
    
    % save results
    G.Ffit_distrib{gi} = Fdist;
    G.details{gi} = details;
      
  end, fprintf('\n'); % next gene

elseif strcmpi(P.method, 'newmethod1_failed_attempt2')

  P=impose_default_value(P,'newmethod1_Frange_to_test','*required*');
  P=impose_default_value(P,'newmethod1_num_neighbors','*required*');
  P=impose_default_value(P,'newmethod1_min_dist2','*required*');

  Frange = P.newmethod1_Frange_to_test;
  mindist2 = P.newmethod1_min_dist2;

  demand_fields(G,{'name','nnon','nsil','Nnon','Nsil'});
  demand_fields(V,'val');

  globalrate_non = sum(G.nnon)/sum(G.Nnon);
  globalrate_sil = sum(G.nsil)/sum(G.Nsil);

 % first convert each covariate to Z-scores and replace missing values with Z=0
  Z = nan(ng,nv);
  for vi=1:nv
    V.missing{vi,1} = (isnan(V.val{vi}) | isinf(V.val{vi}));
    mn = mean(V.val{vi}(~V.missing{vi}));
    sd = std(V.val{vi}(~V.missing{vi}));
    Z(:,vi) = (V.val{vi}-mn)./sd;
    Z(V.missing{vi},vi) = 0;
  end

  % now, for each gene, find distribution on F
  G.Ffit_distribs = cell(ng,1);
  G.Ffit_distrib = cell(ng,1);
  G.Ffit_Q = cell(ng,1);

  for gi=1:ng, if ~mod(gi,1000), fprintf('%d/%d ',gi,ng); end
    if ~analyze_this_gene(gi), continue; end

    % compute and sort squared-distances
    dist2 = sum(bsxfun(@minus,Z,Z(gi,:)).^2,2);
    [tmp ord] = sort(dist2);
    ord(ord==gi)=[]; ord = [gi;ord]; % make sure the gene-itself is listed first
    dist2 = max(dist2,mindist2);

    % for each neighborhood (gene_only  --> up to -->  num_neighbors)
    % compute two distributions:
    %      based on silent
    %      based on nonsilent (excluding the gene itself)
    % for each distribution,
    %      compute the Q = quality = sum(N).^2 / max(Distance2) over genes

    nmuttypes=2;
    Fdist = cell(P.newmethod1_num_neighbors+1,nmuttypes);
    Q = nan(P.newmethod1_num_neighbors+1,nmuttypes);

    for ni=1:P.newmethod1_num_neighbors+1
      for muttype=1:nmuttypes
        Fdist{ni,muttype} = zeros(size(Frange));
        if muttype==1
          type = 'silent';
          gidx = ord(1:ni);
          n = G.nsil(gidx);
          N = G.Nsil(gidx);
          globalrate = globalrate_sil;

          gidx = ord(2:ni);
          type = 'nonsilent';
          n = G.nnon(gidx);
          N = G.Nnon(gidx);
          globalrate = globalrate_non;
        end
        sN=sum(N); sn=sum(n);
        if sN>sn
          lk = binopdf(sn,sN,globalrate*Frange);
          Fdist{ni,muttype} = lk/sum(lk);
          Q(ni,muttype) = sum(N)/(max(dist2(gidx)).^2);
        end
      end % next mutation type
    end % next neighborhood

    % save results
    G.Ffit_distribs{gi,1} = Fdist;
    G.Ffit_Q{gi,1} = Q;
    [tmp idx] = max(Q(:));
    G.Ffit_distrib{gi,1} = Fdist{idx};

  end, fprintf('\n'); % next gene


elseif strcmpi(P.method, 'newmethod1_failed_attempt3')

  P=impose_default_value(P,'newmethod1_num_neighbors','*required*');
  P=impose_default_value(P,'newmethod1_min_dist2','*required*');
  P=impose_default_value(P,'newmethod1_Frange_to_test','*required*');

  mindist2 = P.newmethod1_min_dist2;
  Frange = P.newmethod1_Frange_to_test;

  demand_fields(G,{'name','nnon','nsil','Nnon','Nsil'});
  demand_fields(V,'val');

  globalrate_non = sum(G.nnon)/sum(G.Nnon);
  globalrate_sil = sum(G.nsil)/sum(G.Nsil);

 % first convert each covariate to Z-scores and replace missing values with Z=0
  Z = nan(ng,nv);
  for vi=1:nv
    V.missing{vi,1} = (isnan(V.val{vi}) | isinf(V.val{vi}));
    mn = mean(V.val{vi}(~V.missing{vi}));
    sd = std(V.val{vi}(~V.missing{vi}));
    Z(:,vi) = (V.val{vi}-mn)./sd;
    Z(V.missing{vi},vi) = 0;
  end

  % now, for each gene, find distribution on F
  G.fitrate = nan(ng,1);
  G.nfit = nan(ng,1);
  G.Nfit = nan(ng,1);
  G.Ffit_distrib = cell(ng,1);
  G.details = cell(ng,1);

  for gi=1:ng, if ~mod(gi,1000), fprintf('%d/%d ',gi,ng); end
    if ~analyze_this_gene(gi), continue; end

    % compute and sort squared-distances
    dist2 = sum(bsxfun(@minus,Z,Z(gi,:)).^2,2);
    [tmp ord] = sort(dist2);
    ord(ord==gi)=[]; ord = [gi;ord]; % make sure the gene-itself is listed first
    dist2 = max(dist2,mindist2);

    %%% method:  compute a single nfit and Nfit as follows:
    %%% for each gene in the neighborhood (size specified exactly)
    %%%    for each mutation type (silent for all, nonsilent for all except the gene itself)
    %%%        weight that gene's n and N by distance, by dividing by (dist2/mindist2)
    %%%        (also scale nsil to nnon by ratio of globalrates)
    %%%        add to the running total

    ntot=0; Ntot=0;
    details = {}; didx = 1;

    for ni=0:P.newmethod1_num_neighbors
      if ni==0
        gidx = gi;
        dist2_g = mindist2;
      else
        gidx = ord(ni);
        dist2_g = dist2(gidx);
      end
      name = G.name{gidx};

      for muttype=1:2

        if muttype==1
          type = 'silent';
          n = G.nsil(gidx);
          n = n * (globalrate_non/globalrate_sil);
          N = G.Nsil(gidx);
          globalrate = globalrate_sil;
        elseif muttype==2
          type = 'nonsilent';
          if ni==0, continue; end  % don't count nonsilent mutations of the gene itself
          n = G.nnon(gidx);
          N = G.Nnon(gidx);
          globalrate = globalrate_non;
        end

        if (N>n)
          factor = mindist2/dist2_g;
          ntot=ntot+n*factor;
          Ntot=Ntot+N*factor;
          details.source{didx,1} = type;
          details.ni(didx,1) = ni;
          details.name{didx,1} = name;
          details.dist2(didx,1) = dist2_g;
          details.factor(didx,1) = factor;
          details.n(didx,1) = n;
          details.N(didx,1) = N;
          details.F(didx,1) = (n./N)/globalrate;
          details.ntot(didx,1) = ntot;
          details.Ntot(didx,1) = Ntot;
          didx=didx+1;
        end
      end % next mutation type
    end % next neighbor

    % save results
    G.nfit(gi) = round(ntot);
    G.Nfit(gi) = round(Ntot);
    G.Ffit_distrib{gi} = binopdf(round(ntot),round(Ntot),globalrate_non*Frange);
    G.details{gi} = details;

  end, fprintf('\n'); % next gene

elseif strcmpi(P.method, 'bounded_neighborhood')





end

if isfield(G,'nobs') && isfield(G,'Nobs') && isfield(G,'nfit') && isfield(G,'Nfit') && isfield(G,'Ffit')
  totrate_obs = nansum(G.nobs)/nansum(G.Nobs);
  G.rate_obs = G.nobs./G.Nobs;
  G.Fobs = G.rate_obs / totrate_obs;
  if P.generate_plot, plot_Ffit_vs_Fobs(G,P); end
  if nargout>1
    relrate_adj = as_column(G.Fobs) ./ G.Ffit;
    x = log2(relrate_adj); x(isnan(x)|isinf(x)) = [];
    m = mean(x); s = std(x); v = s^2;
  end
else
  if P.generate_plot
    fprintf('need G.nobs and G.Nobs and G.nfit and G.Nfit in order to generate plot\n');
  end
  if nargout>1
    fprintf('need G.nobs and G.Nobs and G.nfit and G.Nfit in order to calculate m and v\n');
    m=nan; v=nan;
  end
end


