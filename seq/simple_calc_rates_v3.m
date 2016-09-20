function M = simple_calc_rates_v3(M,P)
% modified to use canonical M structure
% modified to use patient subset

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'impute_full_coverage',false);

ng = slength(M.gene);
ncat = M.TOT;
np = slength(M.patient);

P = impose_default_value(P,'bkgd_pat_subset',1:np);
P = impose_default_value(P,'signal_pat_subset',1:np);
bidx = P.bkgd_pat_subset;
sidx = P.signal_pat_subset;
if ischar(bidx), bidx = {bidx}; end
if iscellstr(bidx), bidx = listmap(bidx,M.patient.name); end
if ~isnumeric(bidx), error('invalid P.bidx'); end
if ischar(sidx), sidx = {sidx}; end
if iscellstr(sidx), sidx = listmap(sidx,M.patient.name); end
if ~isnumeric(sidx), error('invalid P.sidx'); end

% if missing coverage, impute full coverage
% for now, use lazy approach of actually building the coverage matrices
% TO DO: switch to intelligent approach of pulling from terr fields on the fly
types = {'sil','non','flank','intron'};
for i=1:length(types)
  covfld = ['N_' types{i} '_cov'];
  terrfld = ['N_' types{i} '_terr'];
  if isfield(M,terrfld) && (~isfield(M,covfld)|| P.impute_full_coverage)
    fprintf('Imputing full coverage for %s',types{i});
    M = setfield(M,covfld,repmat(getfield(M,terrfld),[1 1 np]));
  end
end

% delete any previous results
M = rmfield(M,grep('_c$',fieldnames(M)));
M.gene = rmfield(M.gene,grep('_c$',fieldnames(M.gene)));
flds = {'pfish','nsignal','Nsignal','nbkgd','Nbkgd','qfish','p','q'};
M.gene = rmfield_if_exist(M.gene,flds);

M.rates_c = '------------------------';

haveflank = isfield(M,'N_flank_cov') && isfield(M,'n_flank');
haveintron = isfield(M,'N_intron_cov') && isfield(M,'n_intron');

% CALCULATE GLOBAL RATES
if haveintron
  if size(M.N_intron_cov,3)>1
    M.Nintron_c = round(nansum(sum(M.N_intron_cov(:,:,bidx),3),1));
  else
    M.Nintron_c = round(nansum(M.N_intron_cov,1));
  end
else
  M.Nintron_c = nan(1,ncat);
end
if haveflank
  M.Nflank_c = round(nansum(nansum(M.N_flank_cov(:,:,bidx),3),1));
else
  M.Nflank_c = nan(1,ncat);
end
M.Nsil_c = round(nansum(nansum(M.N_sil_cov(:,:,bidx),3),1));
M.Nsf_c = M.Nflank_c + M.Nsil_c;
M.Nnon_c = round(nansum(nansum(M.N_non_cov(:,:,sidx),3),1));

if haveintron
  if size(M.n_intron,3)>1
    M.nintron_c = round(nansum(nansum(M.n_intron(:,:,bidx),3),1));
  else
    M.nintron_c = round(nansum(M.n_intron,1));
  end
else
  M.nintron_c = nan(1,ncat);
end
if haveflank
  M.nflank_c = round(nansum(nansum(M.n_flank(:,:,bidx),3),1));
else
  M.nflank_c = nan(1,ncat);
end
M.nsil_c = round(nansum(nansum(M.n_silent(:,:,bidx),3),1));
M.nsf_c = M.nflank_c + M.nsil_c;
M.nnon2_c = round(nansum(nansum(M.n_nonsilent_ignoring_null_categ(:,:,sidx),3),1));
M.nnon_c = round(nansum(nansum(M.n_nonsilent(:,:,sidx),3),1));



% rates
M.rintron_c = M.nintron_c./M.Nintron_c;
M.rflank_c = M.nflank_c./M.Nflank_c;
M.rsil_c = M.nsil_c./M.Nsil_c;
M.rsf_c = M.nsf_c./M.Nsf_c;
M.rnon2_c = M.nnon2_c./M.Nnon_c;
M.rnon_c = M.nnon_c./M.Nnon_c;

% FIND BAGELS
if ~isfield(M.gene,'bagel')   % don't recompute if field already exists
  M = calculate_bagels(M,P)
end

% CALCULATE PER-GENE RATES
M.gene.Nbagelflank_c = zeros(ng,ncat);   % (filled in below)
M.gene.Nbagelsil_c = zeros(ng,ncat);
if haveintron
  if size(M.N_intron_cov,3)>1
    M.gene.Nintron_c = squeeze(round(nansum(M.N_intron_cov(:,:,bidx),3)));
  else
    M.gene.Nintron_c = squeeze(round(M.N_intron_cov));
  end
else
  M.gene.Nintron_c = nan(ng,ncat);
end
if haveflank
  M.gene.Nflank_c = squeeze(round(nansum(M.N_flank_cov(:,:,bidx),3)));
else
  M.gene.Nflank_c = nan(ng,ncat);
end
M.gene.Nsil_c = squeeze(round(nansum(M.N_sil_cov(:,:,bidx),3)));
M.gene.Nsf_c = M.gene.Nflank_c + M.gene.Nsil_c;
M.gene.Nsphere_c = nan(ng,ncat);
M.gene.Nnon_c = squeeze(round(nansum(M.N_non_cov(:,:,sidx),3)));

M.gene.n_c = repmat({'---------------------'},ng,1);

M.gene.nbagelflank_c = zeros(ng,ncat);
M.gene.nbagelsil_c = zeros(ng,ncat);
if haveintron
  if size(M.n_intron,3)>1
    M.gene.nintron_c = squeeze(round(nansum(M.n_intron(:,:,bidx),3)));
  else
    M.gene.nintron_c = squeeze(round(M.n_intron));
  end
else
  M.gene.nintron_c = nan(ng,ncat);
end
if haveflank
  M.gene.nflank_c = squeeze(round(nansum(M.n_flank(:,:,bidx),3)));
else
  M.gene.nflank_c = nan(ng,ncat);
end
M.gene.nsil_c = squeeze(round(nansum(M.n_silent(:,:,bidx),3)));
M.gene.nsf_c = M.gene.nflank_c + M.gene.nsil_c;
M.gene.nsphere_c = nan(ng,ncat);
M.gene.nnon_c = squeeze(round(nansum(M.n_nonsilent(:,:,sidx),3)));

for g=1:ng   % fill in information about bagels
  bagel = M.gene.bagel(g,:);
  bagel(isnan(bagel)) = [];
  if isempty(bagel), continue; end
  M.gene.Nbagelflank_c(g,:) = nansum(M.gene.Nflank_c(bagel,:),1);
  M.gene.Nbagelsil_c(g,:) = nansum(M.gene.Nsil_c(bagel,:),1);
  M.gene.nbagelflank_c(g,:) = nansum(M.gene.nflank_c(bagel,:),1);
  M.gene.nbagelsil_c(g,:) = nansum(M.gene.nsil_c(bagel,:),1);
end

% sphere totals
M.gene.Nsphere_c = M.gene.Nbagelflank_c + M.gene.Nbagelsil_c + M.gene.Nflank_c + M.gene.Nsil_c;
M.gene.nsphere_c = M.gene.nbagelflank_c + M.gene.nbagelsil_c + M.gene.nflank_c + M.gene.nsil_c;

M.gene.f_c = repmat({'---------------------'},ng,1);

% CALCULATE CONSERVATIVE EFFECT SIZES

method = 1;

if method==4
  fs = shiftdim(geoseries(0.1,1000,100),-2);
  nfs = length(fs);

  n_sgnl = M.gene.nnon_c;
  N_sgnl = M.gene.Nnon_c;
  ro_sgnl = repmat(M.rnon_c,ng,1);
  N_sgnl(N_sgnl==0) = nan;
  ro_sgnl(ro_sgnl==0) = nan;

  n_bkgd = cat(2,M.gene.nbagelflank_c,M.gene.nbagelsil_c,M.gene.nflank_c,M.gene.nsil_c,M.gene.nsf_c,M.gene.nsphere_c);
  N_bkgd = cat(2,M.gene.Nbagelflank_c,M.gene.Nbagelsil_c,M.gene.Nflank_c,M.gene.Nsil_c,M.gene.Nsf_c,M.gene.Nsphere_c);
  ro_bkgd = repmat(cat(2,M.rflank_c,M.rsil_c,M.rflank_c,M.rsil_c,M.rsf_c,M.rsf_c),ng,1);
  N_bkgd(N_bkgd==0) = nan;
  ro_bkgd(ro_bkgd==0) = nan;


  r_sgnl = bsxfun(@times,ro_sgnl,fs);
  r_bkgd = bsxfun(@times,ro_bkgd,fs);

  p_sgnl = binopdf(repmat(n_sgnl,[1 1 nfs]),repmat(N_sgnl,[1 1 nfs]),r_sgnl);
  p_sgnl_n = bsxfun(@rdivide,p_sgnl,sum(p_sgnl,3));

  p_bkgd = binopdf(repmat(n_bkgd,[1 1 nfs]),repmat(N_bkgd,[1 1 nfs]),r_bkgd);
  p_bkgd_n = bsxfun(@rdivide,p_bkgd,sum(p_bkgd,3));

  submethod = 1;

  if submethod==1

    g=find(strcmp('PTEN',M.gene.name));
    q_bkgd = p_bkgd_n(g,:,:); q_sgnl = p_sgnl_n(g,:,:);
    z_bkgd = squeeze(nansum(q_bkgd,2)); z_sgnl = squeeze(nansum(q_sgnl,2));
    effect = bsxfun(@rdivide,squeeze(fs),squeeze(fs)');
    prob = bsxfun(@times,z_sgnl,z_bkgd');
    xx = zeros(length(fs),1);
    for i=1:length(fs)
      if i==length(fs), idx=find(effect(:)>=fs(end)); else idx=find(effect(:)>=fs(i) & effect(:)<fs(i+1)); end
      xx(i) = sum(prob(idx));
    end
    clf,subplot(4,1,1), hold on, plot(squeeze(q_bkgd)');
    plot(z_bkgd/max(z_bkgd)*max(q_bkgd(:)),'linewidth',2,'color',[0 0 0]); hold off
    set(gca,'xtick',1:5:length(fs),'xticklabel',format_number(squeeze(fs(1:5:end)),3,3)); xlabelvert;
    subplot(4,1,2), hold on, plot(squeeze(q_sgnl)');
    plot(z_sgnl/max(z_sgnl)*max(q_sgnl(:)),'linewidth',2,'color',[1 0 0]); hold off
    set(gca,'xtick',1:5:length(fs),'xticklabel',format_number(squeeze(fs(1:5:end)),3,3)); xlabelvert;
    subplot(4,1,3);hold on
    plot(z_bkgd/max(z_bkgd)/1.1,'color',[0 0 0],'linewidth',2);
    plot(z_sgnl/max(z_sgnl)/1.1,'color',[1 0 0],'linewidth',2);    
    set(gca,'xtick',1:5:length(fs),'xticklabel',format_number(squeeze(fs(1:5:end)),3,3)); xlabelvert;
    subplot(4,1,4);hold on
    plot(xx,'color',[0 0 1],'linewidth',2);
    set(gca,'xtick',1:5:length(fs),'xticklabel',format_number(squeeze(fs(1:5:end)),3,3)); xlabelvert;
    
  elseif submethod==2
    
    g=find(strcmp('PTEN',M.gene.name));
    q_bkgd = squeeze(p_bkgd_n(g,:,:)); q_sgnl = squeeze(p_sgnl_n(g,:,:));
    
    effect = bsxfun(@rdivide,squeeze(fs),squeeze(fs)');
    xx = zeros(size(q_bkgd,1),length(fs));
    xxi = cell(length(fs),1);
    for i=1:length(fs)
      if i==length(fs), xxi{i}=find(effect(:)>=fs(end));
      else xxi{i}=find(effect(:)>=fs(i) & effect(:)<fs(i+1)); end
    end
    
    % find strongest signal channel
    for si=1:size(q_sgnl,1)
      % (using its strongest background channel)
      for bi=1:size(q_bkgd,1)
        prob = bsxfun(@times,q_sgnl(si,:),q_bkgd(bi,:)');
        for i=1:length(fs), xx(bi,i) = sum(prob(xxi{i})); end
      end
    end,end
    


      xx = zeros(length(fs),1);
      for i=1:length(fs)
        if i==length(fs), idx=find(effect(:)>=fs(end));
        else idx=find(effect(:)>=fs(i) & effect(:)<fs(i+1)); end
        xx(i) = sum(prob(idx));
      end

  








  clf,subplot(4,1,1), hold on, plot(squeeze(q_bkgd)');
  plot(z_bkgd/max(z_bkgd)*max(q_bkgd(:)),'linewidth',2,'color',[0 0 0]); hold off
  set(gca,'xtick',1:5:length(fs),'xticklabel',format_number(squeeze(fs(1:5:end)),3,3)); xlabelvert;
  subplot(4,1,2), hold on, plot(squeeze(q_sgnl)');
  plot(z_sgnl/max(z_sgnl)*max(q_sgnl(:)),'linewidth',2,'color',[1 0 0]); hold off
  set(gca,'xtick',1:5:length(fs),'xticklabel',format_number(squeeze(fs(1:5:end)),3,3)); xlabelvert;
  subplot(4,1,3);hold on
  plot(z_bkgd/max(z_bkgd)/1.1,'color',[0 0 0],'linewidth',2);
  plot(z_sgnl/max(z_sgnl)/1.1,'color',[1 0 0],'linewidth',2);
  set(gca,'xtick',1:5:length(fs),'xticklabel',format_number(squeeze(fs(1:5:end)),3,3)); xlabelvert;
  subplot(4,1,4);hold on
  plot(xx,'color',[0 0 1],'linewidth',2);
  set(gca,'xtick',1:5:length(fs),'xticklabel',format_number(squeeze(fs(1:5:end)),3,3)); xlabelvert;




elseif method==1

  [rbagelflank_c sdbagelflank_c] = ratio_and_sd(M.gene.nbagelflank_c,M.gene.Nbagelflank_c);
  [rbagelsil_c sdbagelsil_c] = ratio_and_sd(M.gene.nbagelsil_c,M.gene.Nbagelsil_c);
  [rintron_c sdintron_c] = ratio_and_sd(M.gene.nintron_c,M.gene.Nintron_c);
  [rflank_c sdflank_c] = ratio_and_sd(M.gene.nflank_c,M.gene.Nflank_c);
  [rsil_c sdsil_c] = ratio_and_sd(M.gene.nsil_c,M.gene.Nsil_c);
  [rsf_c sdsf_c] = ratio_and_sd(M.gene.nsf_c,M.gene.Nsf_c);
  [rsphere_c sdsphere_c] = ratio_and_sd(M.gene.nsphere_c,M.gene.Nsphere_c);
  [rnon_c sdnon_c] = ratio_and_sd(M.gene.nnon_c,M.gene.Nnon_c);

  k = 1.96;
  lo_rbagelflank_c = rbagelflank_c-k*sdbagelflank_c;
  lo_rbagelsil_c = rbagelsil_c-k*sdbagelsil_c;
  lo_rintron_c = rintron_c-k*sdintron_c;
  lo_rflank_c = rflank_c-k*sdflank_c;
  lo_rsil_c = rsil_c-k*sdsil_c;
  lo_rsf_c = rsf_c-k*sdsf_c;
  lo_rsphere_c = rsphere_c-k*sdsphere_c;
  lo_rnon_c = rnon_c-k*sdnon_c;

  fbagelflank_c = max(0,bsxfun(@rdivide,lo_rbagelflank_c,M.rflank_c));
  fbagelsil_c = max(0,bsxfun(@rdivide,lo_rbagelsil_c,M.rsil_c));
  fintron_c = max(0,bsxfun(@rdivide,lo_rintron_c,M.rintron_c));
  fflank_c = max(0,bsxfun(@rdivide,lo_rflank_c,M.rflank_c));
  fsil_c = max(0,bsxfun(@rdivide,lo_rsil_c,M.rsil_c));
  fsf_c = max(0,bsxfun(@rdivide,lo_rsf_c,M.rsf_c));
  fsphere_c = max(0,bsxfun(@rdivide,lo_rsphere_c,M.rsf_c));
  fnon_c = max(0,bsxfun(@rdivide,lo_rnon_c,M.rnon_c));
  
elseif method==1.5
  [rbagelflank_c sdbagelflank_c] = ratio_and_sd(M.gene.nbagelflank_c,M.gene.Nbagelflank_c);
  [rbagelsil_c sdbagelsil_c] = ratio_and_sd(M.gene.nbagelsil_c,M.gene.Nbagelsil_c);
  [rflank_c sdflank_c] = ratio_and_sd(M.gene.nflank_c,M.gene.Nflank_c);
  [rsil_c sdsil_c] = ratio_and_sd(M.gene.nsil_c,M.gene.Nsil_c);
  [rsf_c sdsf_c] = ratio_and_sd(M.gene.nsf_c,M.gene.Nsf_c);
  [rsphere_c sdsphere_c] = ratio_and_sd(M.gene.nsphere_c,M.gene.Nsphere_c);
  [rnon_c sdnon_c] = ratio_and_sd(M.gene.nnon_c,M.gene.Nnon_c);

  k = 1.5;
  lo_rbagelflank_c = rbagelflank_c-k*sdbagelflank_c;
  lo_rbagelsil_c = rbagelsil_c-k*sdbagelsil_c;
  lo_rflank_c = rflank_c-k*sdflank_c;
  lo_rsil_c = rsil_c-k*sdsil_c;
  lo_rsf_c = rsf_c-k*sdsf_c;
  lo_rsphere_c = rsphere_c-k*sdsphere_c;
  lo_rnon_c = rnon_c-k*sdnon_c;

  hi_rbagelflank_c = rbagelflank_c+k*sdbagelflank_c;
  hi_rbagelsil_c = rbagelsil_c+k*sdbagelsil_c;
  hi_rflank_c = rflank_c+k*sdflank_c;
  hi_rsil_c = rsil_c+k*sdsil_c;
  hi_rsf_c = rsf_c+k*sdsf_c;
  hi_rsphere_c = rsphere_c+k*sdsphere_c;
  hi_rnon_c = rnon_c+k*sdnon_c;

  lo_fbagelflank_c = bsxfun(@rdivide,lo_rbagelflank_c,M.rflank_c);
  lo_fbagelsil_c = bsxfun(@rdivide,lo_rbagelsil_c,M.rsil_c);
  lo_fflank_c = bsxfun(@rdivide,lo_rflank_c,M.rflank_c);
  lo_fsil_c = bsxfun(@rdivide,lo_rsil_c,M.rsil_c);
  lo_fsf_c = bsxfun(@rdivide,lo_rsf_c,M.rsf_c);
  lo_fsphere_c = bsxfun(@rdivide,lo_rsphere_c,M.rsf_c);
  lo_fnon_c = bsxfun(@rdivide,lo_rnon_c,M.rnon_c);

  hi_fbagelflank_c = bsxfun(@rdivide,hi_rbagelflank_c,M.rflank_c);
  hi_fbagelsil_c = bsxfun(@rdivide,hi_rbagelsil_c,M.rsil_c);
  hi_fflank_c = bsxfun(@rdivide,hi_rflank_c,M.rflank_c);
  hi_fsil_c = bsxfun(@rdivide,hi_rsil_c,M.rsil_c);
  hi_fsf_c = bsxfun(@rdivide,hi_rsf_c,M.rsf_c);
  hi_fsphere_c = bsxfun(@rdivide,hi_rsphere_c,M.rsf_c);
  hi_fnon_c = bsxfun(@rdivide,hi_rnon_c,M.rnon_c);

  % take conservative edge of c.i. (or 1, if c.i. includes 1)

  fbagelflank_c = nan(ng,ncat);
  fbagelsil_c = nan(ng,ncat);
  fflank_c = nan(ng,ncat);
  fsil_c = nan(ng,ncat);
  fsf_c = nan(ng,ncat);
  fsphere_c = nan(ng,ncat);
  fnon_c = nan(ng,ncat);

  fbagelflank_c(lo_fbagelflank_c>1) = lo_fbagelflank_c(lo_fbagelflank_c>1);
  fbagelsil_c(lo_fbagelsil_c>1) = lo_fbagelsil_c(lo_fbagelsil_c>1);
  fflank_c(lo_fflank_c>1) = lo_fflank_c(lo_fflank_c>1);
  fsil_c(lo_fsil_c>1) = lo_fsil_c(lo_fsil_c>1);
  fsf_c(lo_fsf_c>1) = lo_fsf_c(lo_fsf_c>1);
  fsphere_c(lo_fsphere_c>1) = lo_fsphere_c(lo_fsphere_c>1);
  fnon_c(lo_fnon_c>1) = lo_fnon_c(lo_fnon_c>1);

  fbagelflank_c(hi_fbagelflank_c<1) = hi_fbagelflank_c(hi_fbagelflank_c<1);
  fbagelsil_c(hi_fbagelsil_c<1) = hi_fbagelsil_c(hi_fbagelsil_c<1);
  fflank_c(hi_fflank_c<1) = hi_fflank_c(hi_fflank_c<1);
  fsil_c(hi_fsil_c<1) = hi_fsil_c(hi_fsil_c<1);
  fsf_c(hi_fsf_c<1) = hi_fsf_c(hi_fsf_c<1);
  fsphere_c(hi_fsphere_c<1) = hi_fsphere_c(hi_fsphere_c<1);
  fnon_c(hi_fnon_c<1) = hi_fnon_c(hi_fnon_c<1);
  
elseif method==2

  alpha = 0.05;
  [tmp lowbagelflank_c tmp] = binofit_2d(M.gene.nbagelflank_c,M.gene.Nbagelflank_c,alpha);
  [tmp lowbagelsil_c tmp] = binofit_2d(M.gene.nbagelsil_c,M.gene.Nbagelsil_c,alpha);
  [tmp lowflank_c tmp] = binofit_2d(M.gene.nflank_c,M.gene.Nflank_c,alpha);
  [tmp lowsil_c tmp] = binofit_2d(M.gene.nsil_c,M.gene.Nsil_c,alpha);
  [tmp lowsf_c tmp] = binofit_2d(M.gene.nsf_c,M.gene.Nsf_c,alpha);
  [tmp lowsphere_c tmp] = binofit_2d(M.gene.nsphere_c,M.gene.Nsphere_c,alpha);
  [tmp lownon_c tmp] = binofit_2d(M.gene.nnon_c,M.gene.Nnon_c,alpha);

  rbagelflank_c = lowbagelflank_c;
  rbagelsil_c = lowbagelsil_c;
  rflank_c = lowflank_c;
  rsil_c = lowsil_c;
  rsf_c = lowsf_c;
  rsphere_c = lowsphere_c;
  rnon_c = lownon_c;

  fbagelflank_c = bsxfun(@rdivide,rbagelflank_c,M.rflank_c);% fbagelflank_c(isnan(fbagelflank_c))=0;
  fbagelsil_c = bsxfun(@rdivide,rbagelsil_c,M.rsil_c);% fbagelsil_c(isnan(fbagelsil_c))=0;
  fflank_c = bsxfun(@rdivide,rflank_c,M.rflank_c);% fflank_c(isnan(fflank_c))=0;
  fsil_c = bsxfun(@rdivide,rsil_c,M.rsil_c);% fsil_c(isnan(fsil_c))=0;
  fsf_c = bsxfun(@rdivide,rsf_c,M.rsf_c);% fsf_c(isnan(fsf_c))=0;
  fsphere_c = bsxfun(@rdivide,rsphere_c,M.rsf_c);% fsphere_c(isnan(fsphere_c))=0;
  fnon_c = bsxfun(@rdivide,rnon_c,M.rnon_c);% fnon_c(isnan(fnon_c))=0;


end

%%%%%% SPECIAL CASE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  genes with RPKM=0  --> assess penalty to rnon_c
if isfield(M.gene,'exclude_because_zero')
  penalty = 10;
  rnon_c(M.gene.exclude_because_zero,:) = rnon_c(M.gene.exclude_because_zero,:) / penalty;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if method==4


else
  M.gene.fbagelflank_c = fbagelflank_c;
  M.gene.fbagelsil_c = fbagelsil_c;
  M.gene.fintron_c = fintron_c;
  M.gene.fflank_c = fflank_c;
  M.gene.fsil_c = fsil_c;
  M.gene.fsf_c = fsf_c;
  M.gene.fsphere_c = fsphere_c;
  M.gene.fnon_c = fnon_c;
end

% CALCULATE FINAL SCORES

M.gene.x_c = repmat({'---------------------'},ng,1);

if method==4

else

  bkgd = [fbagelflank_c fbagelsil_c fintron_c fflank_c fsil_c fsf_c fsphere_c];
  signal = [fnon_c];

  [bkgd bkgd_idx] = max(bkgd,[],2);
  [signal signal_idx] = max(signal,[],2);

  bkgd(bkgd<0.1) = 0.1;
  bkgd(isnan(bkgd)) = 1;
  bkgd(isinf(bkgd)) = 10;

  signal(isnan(signal)) = 0;
  signal(isinf(signal)) = 10;

  M.gene.bkgd_c = bkgd;
  M.gene.signal_c = signal;

  score_c = (signal./bkgd);

  % Fisher p-value

  nbkgd = [M.gene.nbagelflank_c M.gene.nbagelsil_c M.gene.nintron_c M.gene.nflank_c M.gene.nsil_c M.gene.nsf_c M.gene.nsphere_c];
  nsignal = [M.gene.nnon_c];

  Nbkgd = [M.gene.Nbagelflank_c M.gene.Nbagelsil_c M.gene.Nintron_c M.gene.Nflank_c M.gene.Nsil_c M.gene.Nsf_c M.gene.Nsphere_c];
  Nsignal = [M.gene.Nnon_c];

  M.gene.pfish = nan(ng,1); M.gene.qfish = nan(ng,1);
  M.gene.nsignal = nan(ng,1); M.gene.Nsignal = nan(ng,1);
  M.gene.nbkgd = nan(ng,1); M.gene.Nbkgd = nan(ng,1);

if 0
  ncol = size(nbkgd,2) / size(fsil_c,2);
  fprintf('Fisher: ');
  for i=1:ng, if ~mod(i,1000), fprintf('%d/%d ',i,ng); end
    col = signal_idx(i);
    sn = nsignal(i,col)*ones(ncol,1);
    sN = Nsignal(i,col)*ones(ncol,1);
    M.gene.nsignal(i) = sn(1);
    M.gene.Nsignal(i) = sN(1);
    bn = nbkgd(i,(col:ncol:end)');
    bN = Nbkgd(i,(col:ncol:end)');
    pf = fisher_exact_test(sn,sN,bn,bN);
    [m idx] = min(pf);
    M.gene.nbkgd(i) = bn(idx);
    M.gene.Nbkgd(i) = bN(idx);
    M.gene.pfish(i) = m;
  end, fprintf('\n');
  M.gene.qfish = calc_fdr_value(M.gene.pfish);
end

end


M.gene.score_c = score_c;

M.gene.y_c = repmat({'---------------------'},ng,1);


%%%%%% SPECIAL CASE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  genes with RPKM=0  --> assess penalty to rnon_c
%if isfield(M.gene,'exclude_because_zero')
%  M.gene.score_c(M.gene.exclude_because_zero,:) = -100;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CALCULATE P-VALUES
% using Empirical Bayes methodology

s = M.gene.score_c;

% s = 0         50-70 %
% 0 < s <= 1    20-40 %
% 1 < s <= 2    2-3 %
% 2 < s <= 3    0.3 - 1 %
% 3 < s <= 4    0.04 - 0.06 %
% 4 < s <= 5    ~5 genes

% s = 0         log10p = 0     p = 1
% s = 1         log10p ~ 0.5   p ~ 0.3
% s = 2         log10p ~ 1     p ~ 0.1
% s = 3         log10p = 4     p ~ 10^-4
% s = 4         log10p ~ 5     p ~ 10^-5
% s>10          log10p = inf   p ~ 0

% functional form:
%    log10p = ks
%    --> fit k to flat Q-Q plot

method=1;

if method==1

  % k=1; lp = s*k; p = 10.^(-lp); 

  % k=1 works just fine for luad337, but need to fit for other cases..
  % TO DO: adjust for number of (p==1) genes

  p = 10.^(-s);

elseif method==2

  s = M.gene.score_c;

  %s(s==0) = min(s(s>0));
  %idx = find(s>0);
  %ls = log2(s(idx));
  %ls = s;
  %ls(isinf(ls)|isnan(ls))=[];
  %[mu sig] = normfit(ls);

  %p = ones(ng,1);
  %p(idx) = 1-normcdf(ls,mu,sig);

  p = 10.^(-s);

  p(p<0) = 0; p(p>1) = 1;
  po=p;
  
  %qq(p); png
  
  a=p;
  y=-log10(a);
  y(isinf(y)) = max(y(~isinf(y)));
  na = length(a);
  x=-log10((1:na)'./(na+1));
  x = sort(x);
  y = sort(y);

  bins = (0:0.1:2.5)';
  ybn = zeros(length(bins)-1,1);
  yb = zeros(length(bins)-1,1);
  for i=1:length(bins)-1
    idx = find(x>=bins(i) & x<bins(i+1));
    ybn(i) = length(idx);
    yb(i) = median(y(idx));
  end
  
  %keyboard
  i1 = find(yb>0,1);
  if isempty(i1), i1=1; end
  zzz = cumsum(ybn)/sum(ybn);
  i2 = find(zzz>0.999,1);  % ignore top 0.1% of genes
  if isempty(i2), i2=length(zzz)+1; end
  ii = i1:(i2-1);

  xfit = bins(ii);
  yfit = yb(ii);
  
  as = -2:0.1:10; na=length(as);
  bs = -2:0.1:4; nb=length(bs);
  r = zeros(na,nb);
  for ai=1:na, for bi=1:nb
      a=as(ai); b=bs(bi);
      yi = a*yfit+b;
      r(ai,bi) = sum(abs(yi-xfit).^2);
  end,end
  [tmp idx] = min(r(:));
  ci = coords([na nb],idx);
  a=as(ci(1)); b=bs(ci(2));
    
  %p = 10.^-(a*-log10(po)+b);
    
  p = 10.^-(a*s+b);

  p(po==1)=1;

end
    
p(p>1) = 1; p(p<0) = 0;

qq(p);ylim([0 20]);
%png

q = calc_fdr_value(p);

M.gene.p = p;
M.gene.q = q;

% reporting structure
M.g = [];
M.g.name = M.gene.name;

M.g.N = M.gene.Nnon_c(:,end);
M.g.nintron = M.gene.nintron_c(:,end);
M.g.nflank = M.gene.nflank_c(:,end);
M.g.nsil = M.gene.nsil_c(:,end);
M.g.nnon = M.gene.nnon_c(:,end);
M.g.nnull = M.gene.nnon_c(:,ncat-1);

nintron = M.gene.nintron_c(:,end); Nintron = M.gene.Nintron_c(:,end);
fintron = (nintron./Nintron)/(sum(nintron)/sum(Nintron)); fintron(Nintron==0) = 0;

nflank = M.gene.nflank_c(:,end); Nflank = M.gene.Nflank_c(:,end);
fflank = (nflank./Nflank)/(sum(nflank)/sum(Nflank)); fflank(Nflank==0) = 0;

nsil = M.gene.nsil_c(:,end); Nsil = M.gene.Nsil_c(:,end);
fsil = (nsil./Nsil)/(sum(nsil)/sum(Nsil)); fsil(Nsil==0) = 0;

nsf = M.gene.nsf_c(:,end); Nsf = M.gene.Nsf_c(:,end);
fsf = (nsf./Nsf)/(sum(nsf)/sum(Nsf)); fsf(Nsf==0) = 0;

nbagelflank = M.gene.nbagelflank_c(:,end); Nbagelflank = M.gene.Nbagelflank_c(:,end);
fbagelflank = (nbagelflank./Nbagelflank)/(sum(nflank)/sum(Nflank)); fbagelflank(Nbagelflank==0) = 0;

nbagelsil = M.gene.nbagelsil_c(:,end); Nbagelsil = M.gene.Nbagelsil_c(:,end);
fbagelsil = (nbagelsil./Nbagelsil)/(sum(nsil)/sum(Nsil)); fbagelsil(Nbagelsil==0) = 0;

nsf = M.gene.nsf_c(:,end); Nsf = M.gene.Nsf_c(:,end);
fsf = (nsf./Nsf)/(sum(nsf)/sum(Nsf)); fsf(Nsf==0) = 0;

nsphere = M.gene.nsphere_c(:,end); Nsphere = M.gene.Nsphere_c(:,end);
fsphere = (nsphere./Nsphere)/(sum(nsf)/sum(Nsf)); fsphere(Nsphere==0) = 0;

nnon = M.gene.nnon_c(:,end); Nnon = M.gene.Nnon_c(:,end);
fnon = (nnon./Nnon)/(sum(nnon)/sum(Nnon)); fnon(Nnon==0) = 0;

M.g.Fbflk = round(fbagelflank);
M.g.Fbsil = round(fbagelsil);
M.g.Fintr = round(fintron);
M.g.Fflk = round(fflank);
M.g.Fsil = round(fsil);
M.g.Fsf = round(fsf);
M.g.Fsphere = round(fsphere);
M.g.Fnon = round(fnon);

%M.g.pfish = M.gene.pfish;
%M.g.qfish = M.gene.qfish;
M.g.score = M.gene.score_c;
M.g.p = M.gene.p;
M.g.q = M.gene.q;

M.g = sort_struct(M.g,{'p','score','nnon'},[1 -1 -1]);





