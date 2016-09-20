function [M,reportstruct] = estimate_mutation_rates(M,P)
% Mike Lawrence 2009-07-02

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'force_recalculate',false);
P=impose_default_value(P,'impute_full_coverage',~isfield(M,'N_cov'));
P=impose_default_value(P,'presumed_driver_genes','');
P=impose_default_value(P,'use_sample_specific_mutation_rates',true);

if isfield(M,'multi')   % multimode
  nsets = length(M.multi);
  reportstruct = cell(nsets+1,1);
  for i=1:nsets
    [M.multi{i} reportstruct{i}] = estimate_mutation_rates(M.multi{i},P);
    reportstruct{i}.isetname = repmat({M.multi{i}.name},slength(reportstruct{i}),1);
  end
  tmp = rmfield(M,'multi');
  [tmp reportstruct{nsets+1}] = estimate_mutation_rates(tmp,P);
  M.mutrate = tmp.mutrate;
  reportstruct{nsets+1}.isetname = repmat({'combined'},slength(reportstruct{nsets+1}),1);
  reportstruct = concat_structs(reportstruct);
  reportstruct = order_fields_first(reportstruct,{'isetname'});
  return
end

if ~P.impute_full_coverage
  ns_s_info_available = (isfield(M,'N_sil_cov') && isfield(M,'N_non_cov'));
else
  fprintf('estimate_mutation_rates: imputing full coverage\n');
  ns_s_info_available = (isfield(M,'N_sil_terr') && isfield(M,'N_non_terr'));
end

reportstruct = [];

if P.force_recalculate || ~isfield(M,'mutrate')

  M.mutrate = [];

  % choose which genes to include
  if ~isempty(P.presumed_driver_genes)
    if iscell(P.presumed_driver_genes)
      P.presumed_driver_genes = concat(P.presumed_driver_genes,'|');
    end
    fprintf('Genes excluded from BMR calculation:\n\t%s\n',...
            regexprep(P.presumed_driver_genes,'|',' '));
  end
  exclude_idx = grep(['^(' P.presumed_driver_genes ')$'],M.gene.name,1);
  disp(M.gene.name(exclude_idx));
  gidx = setdiff(1:M.ng,exclude_idx);

  % calculate global BMR (even if using ssBMR's)
  n = []; N = [];
  for c=1:M.TOT
    if ~P.impute_full_coverage
      N(c,1) = fullsum(M.N_cov(gidx,c,:));
    else
      N(c,1) = fullsum(M.N_terr(gidx,c)) .* M.np;
    end
    n(c,1) = fullsum(M.n_nonsilent(gidx,c,:));
  end
  rate = n./N;
  M.mutrate.n = n(1:M.TOT-1)';
  M.mutrate.N = N(1:M.TOT-1)';
  M.mutrate.rate = rate(1:M.TOT-1)';
  M.mutrate.rel = (rate(1:M.TOT-1)/rate(M.TOT))';
  M.mutrate.tot.hat = rate(M.TOT);

  if ns_s_info_available
    N_sil = []; N_non = [];
    for c=1:M.TOT
      if ~P.impute_full_coverage
        N_sil(c,1) = fullsum(M.N_sil_cov(gidx,c,:));
        N_non(c,1) = fullsum(M.N_non_cov(gidx,c,:));
      else
        N_sil(c,1) = fullsum(M.N_sil_terr(gidx,c)) .* M.np;
        N_non(c,1) = fullsum(M.N_non_terr(gidx,c)) .* M.np;
      end
    end
    exp_ns_s_ratio = N_non./N_sil;
    exp_ns_s_ratio(M.TOT-M.NUM_INDEL_CLASSES:M.TOT-1) = nan;
    M.mutrate.exp_ns_s_ratio = exp_ns_s_ratio(1:M.TOT-1)';
    M.mutrate.tot.exp_ns_s_ratio = exp_ns_s_ratio(M.TOT);

    % also calculate nonsilent rates expected from silent
    n2 = []; N2 = [];
    for c=1:M.TOT
      if ~P.impute_full_coverage
        N2(c,1) = fullsum(M.N_cov(gidx,c,:));
      else
        N2(c,1) = fullsum(M.N_terr(gidx,c)) .* M.np;
      end
      n2(c,1) = fullsum(M.n_silent(gidx,c,:));
    end
    rate2 = n2./N2;
    extrapolated_ns_rate = rate2 .* exp_ns_s_ratio;
    M.mutrate.rate_extrapolated_from_silent = extrapolated_ns_rate';

  end

  if isfield(M,'note1') && strcmp(M.note1,'FOR SILENT P-VALUE CALC')
    tt = 'silent';
  else 
    tt = 'nonsilent';
  end
  fprintf('\nEstimated %s background mutation rates\n\n',tt);
  for c=1:M.TOT
    fprintf('%19s  %9.2d  %5.1f/Mb  (relative = %0.2f)',M.mutclass{c},rate(c),rate(c)*1e6,rate(c)/rate(end));
    if ns_s_info_available
      if ~isnan(exp_ns_s_ratio(c))
        fprintf('  expected NS/S ratio = %3.1f',exp_ns_s_ratio(c));
      end
    end
    fprintf('\n');
  end
  fprintf('\n');
  
  X=[];
  X.category = as_column(M.mutclass);
  X.n = n;
  X.N = N;
  X.rate = rate;
  X.rate_per_mb = rate*1e6;
  X.relative_rate = rate/rate(end);
  if ns_s_info_available
    X.exp_ns_s_ratio = exp_ns_s_ratio;
  end
  reportstruct = X;
    
  if P.use_sample_specific_mutation_rates
    n = []; N = [];
    for p=1:M.np, for c=1:M.TOT
      if ~P.impute_full_coverage
        N(c,p) = fullsum(M.N_cov(gidx,c,p));
      else
        N(c,p) = fullsum(M.N_terr(gidx,c));
      end
      n(c,p) = fullsum(M.n_nonsilent(gidx,c,p));
    end, end
    rate = n./N;
      
    M.mutrate.ss = [];
    M.mutrate.ss.n = n(1:M.TOT-1,:)';
    M.mutrate.ss.N = N(1:M.TOT-1,:)';
    M.mutrate.ss.rate = rate(1:M.TOT-1,:)';
    M.mutrate.ss.rel = bsxfun(@rdivide,rate(1:M.TOT-1,:),rate(M.TOT,:))';
    M.mutrate.ss.tot.hat = rate(M.TOT,:);

    if ns_s_info_available
      N_sil = []; N_non = [];
      for p=1:M.np, for c=1:M.TOT
        if ~P.impute_full_coverage
          N_sil(c,p) = fullsum(M.N_sil_cov(gidx,c,p));
          N_non(c,p) = fullsum(M.N_non_cov(gidx,c,p));
        else
          N_sil(c,p) = fullsum(M.N_sil_terr(gidx,c));
          N_non(c,p) = fullsum(M.N_non_terr(gidx,c));
        end
      end, end
      exp_ns_s_ratio = N_non./N_sil;
      exp_ns_s_ratio(M.TOT-M.NUM_INDEL_CLASSES:M.TOT-1,:) = nan;
      M.mutrate.ss.exp_ns_s_ratio = exp_ns_s_ratio(1:M.TOT-1,:)';
      M.mutrate.ss.tot.exp_ns_s_ratio = exp_ns_s_ratio(M.TOT,:);

      % also calculate nonsilent rates expected from silent
      n2 = []; N2 = [];
      for p=1:M.np, for c=1:M.TOT
        if ~P.impute_full_coverage
          N2(c,p) = fullsum(M.N_cov(gidx,c,p));
        else
          N2(c,p) = fullsum(M.N_terr(gidx,c));
        end
        n2(c,p) = fullsum(M.n_silent(gidx,c,p));
      end,end
      rate2 = n2./N2;
      extrapolated_ns_rate = rate2 .* exp_ns_s_ratio;
      M.mutrate.ss.rate_extrapolated_from_silent = extrapolated_ns_rate(1:M.TOT-1,:)';
      M.mutrate.ss.tot.hat = extrapolated_ns_rate(M.TOT,:);
    end
     
    fprintf('\nEstimated SAMPLE-SPECIFIC nonsilent background mutation rates\n\n');
    for p=1:M.np
      fprintf('%s\n',M.patient.name{p});
      for c=1:M.TOT
        fprintf('%19s  %9.2d  %5.1f/Mb  (relative = %0.2f)',M.mutclass{c},rate(c,p),rate(c,p)*1e6,rate(c,p)/rate(end,p));
        if ns_s_info_available
          if ~isnan(exp_ns_s_ratio(c,p))
            fprintf('  expected NS/S ratio = %3.1f',exp_ns_s_ratio(c,p));
          end
        end
        fprintf('\n');
      end
      fprintf('\n');
    end

%    X=[];
%    tmp = repmat(M.patient.name',M.TOT,1);
%    X.patient = tmp(:);
%    X.category = repmat(as_column(M.mutclass),M.np,1);
%    X.n = n(:);
%    X.N = N(:);
%    X.rate = rate(:);
%    X.rate_per_mb = rate(:)*1e6;
%    tmp = bsxfun(@rdivide,rate,rate(M.TOT,:));
%    X.relative_rate = tmp(:);
%    if ns_s_info_available
%      X.exp_ns_s_ratio = exp_ns_s_ratio(:);
%    end   
%    reportstruct = X;
  end
  
else
  fprintf('estimate_mutation_rates: using analysis already performed\n');
end
