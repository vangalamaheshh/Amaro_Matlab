function n_work = estimate_nonsilent_passengers_from_silent(n_work,M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'impute_full_coverage',~isfield(M,'N_sil_cov'));
P=impose_default_value(P,'pick_random_passengers',false);

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

fprintf('estimate_nonsilent_passengers_from_silent\n');
silent = true;

demand_fields(M,{'mutrate','n_nonsilent_ignoring_null_categ','n_silent'});
if ~P.impute_full_coverage
  demand_fields(M,{'N_sil_cov','N_non_cov'});
else
  demand_fields(M,{'N_sil_terr','N_non_terr'});
end

ncat1 = M.TOT - M.NUM_INDEL_CLASSES - 1;
ncat2 = M.TOT - 1;

for g=gta

  gname = M.gene.name{g};

  ncat = ncat1;       % if <10 silent mutations, then don't allow removal of null mutations
  nsil_obs = sum(sum(M.n_silent(g,1:ncat,:),3));
  nnon_obs = sum(sum(M.n_nonsilent_ignoring_null_categ(g,1:ncat,:),3));
  nnon_obs2 = sum(sum(M.n_nonsilent(g,1:ncat,:),3));
  if (nsil_obs>=10)   % if >=10 silent mutations, then allow removal of null mutations
    ncat = ncat2;
    nsil_obs = sum(sum(M.n_silent(g,1:ncat,:),3));
    nnon_obs = sum(sum(M.n_nonsilent_ignoring_null_categ(g,1:ncat,:),3));
    nnon_obs2 = sum(sum(M.n_nonsilent(g,1:ncat,:),3));
  end

  if nsil_obs<2 || nnon_obs2==0, continue; end

  ns_s_ratio_obs = nnon_obs / nsil_obs;

  % calculate expected non/sil ratio

%  Nnon = sum(M.N_non_cov(g,M.TOT,:),3);
%  Nsil = sum(M.N_sil_cov(g,M.TOT,:),3);
%  ns_s_ratio_exp = Nnon / Nsil;

   if ~P.impute_full_coverage
     Nnon = sum(M.N_non_cov(g,1:ncat,:),3);
     Nsil = sum(M.N_sil_cov(g,1:ncat,:),3);
   else
     Nnon = M.N_non_terr(g,1:ncat).*M.np;
     Nsil = M.N_sil_terr(g,1:ncat).*M.np;
   end
   ns_s_ratio_exp = weighted_mean(Nnon./Nsil,M.mutrate.rel(1:ncat));

  % given observed number of silent mutations and fnon_exp, calculate expected number of nonsilent
  if nsil_obs<5
    nnon_exp = (nsil_obs-1) * ns_s_ratio_exp;     % -1 to be generous
  elseif nsil_obs<10
    nnon_exp = nsil_obs * ns_s_ratio_exp;
  elseif nsil_obs<30
    nnon_exp = nsil_obs * 1.2 * ns_s_ratio_exp;
  else
    nnon_exp = nsil_obs * 1.5 * ns_s_ratio_exp;
  end
  npass = floor(nnon_exp);

  % nominate passengers and remove

  if P.pick_random_passengers   % (actually less effective at removing bad genes e.g. TTN, in projection mode)
    cidx = randperm(ncat);
    pidx = randperm(size(n_work,3));
  else
    cidx = 1:ncat;
    pidx = 1:size(n_work,3);
  end

  nremoved = 0;
  while(npass>0)
    [c p] = find(squeeze(n_work(g,cidx,pidx))>0,1);
    if isempty(p), break; end
    n_work(g,cidx(c),pidx(p)) = n_work(g,cidx(c),pidx(p))-1;
    npass=npass-1;
    nremoved=nremoved+1;
  end

  if ~silent
    if nremoved>0
      fprintf('\tgene %-10s \t(NS/S exp %0.2f  obs %d/%d = %0.2f)\t-->  %d  nonsilent passengers removed\n',...
              gname,ns_s_ratio_exp,nnon_obs,nsil_obs,ns_s_ratio_obs,nremoved);
    end
  end

end




