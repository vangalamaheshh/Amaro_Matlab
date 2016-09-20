function M = simple_calc_rates(M,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'vidx', [3 6 10 17 15]);
P = impose_default_value(P,'bagel_size', 30);

% if doesn't exist already, create a "total" category
if ~strcmpi(M.categ.name(end),'Total')
  fprintf('Adding "total" category\n');
  totcateg = reorder_struct(M.categ,'end');
  totcateg.name = {'Total'}; totcateg.type = {'Total'};
  M.categ = concat_structs_keep_all_fields({M.categ,totcateg});
  M.gene_terr(:,end+1) = M.gene_terr(:,end);
  M.gene_sil_terr(:,end+1) = M.gene_sil_terr(:,end);
  M.gene_non_terr(:,end+1) = M.gene_non_terr(:,end);
  M.gene_flank_terr(:,end+1) = M.gene_flank_terr(:,end);
  M.gene_cov(:,:,end+1) = M.gene_cov(:,:,end);
  M.gene_sil_cov(:,:,end+1) = M.gene_sil_cov(:,:,end);
  M.gene_non_cov(:,:,end+1) = M.gene_non_cov(:,:,end);
  M.gene_flank_cov(:,:,end+1) = M.gene_flank_cov(:,:,end);
  M.ncat = slength(M.categ);
end

ng = slength(M.gene);
ncat = slength(M.categ);
np = slength(M.sample);

M = rmfield(M,grep('_c$',fieldnames(M)));
M.gene = rmfield(M.gene,grep('_c$',fieldnames(M.gene)));

% total
M.Nflank_c = squeeze(round(sum(sum(M.gene_flank_cov,1),2)))';
M.Nsil_c = squeeze(round(sum(sum(M.gene_sil_cov,1),2)))';
M.Nsf_c = M.Nflank_c + M.Nsil_c;
M.Nnon_c = squeeze(round(sum(sum(M.gene_non_cov,1),2)))';
idx = find(M.mut.is_flank); M.nflank_c = histc(M.mut.categ(idx),1:ncat)';
idx = find(M.mut.is_coding & M.mut.is_silent); M.nsil_c = histc(M.mut.categ(idx),1:ncat)';
idx = find(M.mut.is_coding & ~M.mut.is_silent);
M.nsf_c = M.nflank_c + M.nsil_c;
M.nnon2_c = histc(M.mut.categ_ignoring_null_categ(idx),1:ncat)';
M.nnon_c = histc(M.mut.categ(idx),1:ncat)';

M.nflank_c(:,end) = sum(M.nflank_c,2);
M.nsil_c(:,end) = sum(M.nsil_c,2);
M.nsf_c(:,end) = sum(M.nsf_c,2);
M.nnon2_c(:,end) = sum(M.nnon2_c,2);
M.nnon_c(:,end) = sum(M.nnon_c,2);

M.rflank_c = M.nflank_c./M.Nflank_c;
M.rsil_c = M.nsil_c./M.Nsil_c;
M.rsf_c = M.nsf_c./M.Nsf_c;
M.rnon2_c = M.nnon2_c./M.Nnon_c;
M.rnon_c = M.nnon_c./M.Nnon_c;

% FIND BAGELS

if ~isfield(M.gene,'bagel')   % don't recompute if field already exists
  demand_field(M,'V');
  nv = slength(M.V);

  % convert to Z-scores
  Z = nan(ng,nv);
  for vi=1:nv
    if length(M.V.val{vi})~=ng, error('V.val must match genelist'); end
    missing = (isnan(M.V.val{vi}) | isinf(M.V.val{vi}));
    mn = mean(M.V.val{vi}(~missing));
    sd = std(M.V.val{vi}(~missing));
    Z(:,vi) = (M.V.val{vi}-mn)./sd;
    Z(missing,vi) = nan;
  end

  nb = P.bagel_size;
  M.gene.bagel = nan(ng,nb);
  fprintf('Finding bagels:');
  for gi=1:ng, if ~mod(gi,1000), fprintf(' %d/%d',gi,ng); end
 
    gz = Z(gi,P.vidx);
    idx = find(~isnan(gz));
    if isempty(idx), continue; end   % no data = no bagel!
    vidx = P.vidx(idx);

    others = [1:(gi-1) (gi+1):ng];
    dist2 = bsxfun(@minus,Z(others,vidx),Z(gi,vidx)).^2;

    % find nearest neighbors
    method = 2;
    if method == 1

      dist2 = nanmean(dist2,2);
      [tmp,ord] = sort(dist2);

    elseif method == 2
      % favor genes with more information
      %  each category can contribute at most 10 points; this is downweighted by dist^2

      pts = 10*ones(size(dist2));
      pts(dist2>0.001) = 9;
      pts(dist2>0.003) = 8;
      pts(dist2>0.01) = 7;
      pts(dist2>0.03) = 6;
      pts(dist2>0.1) = 5;
      pts(dist2>0.3) = 4;
      pts(dist2>1) = 3;
      pts(dist2>2) = 2;
      pts(dist2>3) = 1;
      pts(dist2>4) = 0;
      pts(isnan(dist2)) = 0;

      pts = sum(pts,2);
      [tmp,ord] = sort(pts,'descend');

    end

    % find nearest neighbors
    bagel = others(ord(1:nb));
    M.gene.bagel(gi,:) = bagel;
  end, fprintf('\n');
  
end

% per gene
M.gene.Nbagelflank_c = zeros(ng,ncat);   % (filled in below)
M.gene.Nbagelsil_c = zeros(ng,ncat);
M.gene.Nflank_c = squeeze(round(sum(M.gene_flank_cov,2)));
M.gene.Nsil_c = squeeze(round(sum(M.gene_sil_cov,2)));
M.gene.Nsf_c = M.gene.Nflank_c + M.gene.Nsil_c;
M.gene.Nsphere_c = nan(ng,ncat);
M.gene.Nnon_c = squeeze(round(sum(M.gene_non_cov,2)));

M.gene.n_c = repmat({'---------------------'},ng,1);

M.gene.nbagelflank_c = zeros(ng,ncat);
M.gene.nbagelsil_c = zeros(ng,ncat);

idx = find(M.mut.is_flank);
M.gene.nflank_c = hist2d_fast(M.mut.gene_idx(idx),M.mut.categ(idx),1,ng,1,ncat);
idx = find(M.mut.is_coding & M.mut.is_silent);
M.gene.nsil_c = hist2d_fast(M.mut.gene_idx(idx),M.mut.categ(idx),1,ng,1,ncat);
idx = find(M.mut.is_coding & ~M.mut.is_silent);
M.gene.nsf_c = M.gene.nflank_c + M.gene.nsil_c;
M.gene.nsphere_c = nan(ng,ncat);
M.gene.nnon_c = hist2d_fast(M.mut.gene_idx(idx),M.mut.categ(idx),1,ng,1,ncat);

M.gene.nflank_c(:,end) = sum(M.gene.nflank_c,2);
M.gene.nsil_c(:,end) = sum(M.gene.nsil_c,2);
M.gene.nsf_c(:,end) = sum(M.gene.nsf_c,2);
M.gene.nnon_c(:,end) = sum(M.gene.nnon_c,2);

for g=1:ng   % fill in information about bagels
  bagel = M.gene.bagel(g,:);
  bagel(isnan(bagel)) = [];
  if isempty(bagel), continue; end
  M.gene.Nbagelflank_c(g,:) = sum(M.gene.Nflank_c(bagel,:),1);
  M.gene.Nbagelsil_c(g,:) = sum(M.gene.Nsil_c(bagel,:),1);
  M.gene.nbagelflank_c(g,:) = sum(M.gene.nflank_c(bagel,:),1);
  M.gene.nbagelsil_c(g,:) = sum(M.gene.nsil_c(bagel,:),1);
end

% sphere totals
M.gene.Nsphere_c = M.gene.Nbagelflank_c + M.gene.Nbagelsil_c + M.gene.Nflank_c + M.gene.Nsil_c;
M.gene.nsphere_c = M.gene.nbagelflank_c + M.gene.nbagelsil_c + M.gene.nflank_c + M.gene.nsil_c;

M.gene.f_c = repmat({'---------------------'},ng,1);

%%%% CONSERVATIVE EFFECT SIZES

[rbagelflank_c sdbagelflank_c] = ratio_and_sd(M.gene.nbagelflank_c,M.gene.Nbagelflank_c);
[rbagelsil_c sdbagelsil_c] = ratio_and_sd(M.gene.nbagelsil_c,M.gene.Nbagelsil_c);
[rflank_c sdflank_c] = ratio_and_sd(M.gene.nflank_c,M.gene.Nflank_c);
[rsil_c sdsil_c] = ratio_and_sd(M.gene.nsil_c,M.gene.Nsil_c);
[rsf_c sdsf_c] = ratio_and_sd(M.gene.nsf_c,M.gene.Nsf_c);
[rsphere_c sdsphere_c] = ratio_and_sd(M.gene.nsphere_c,M.gene.Nsphere_c);
[rnon_c sdnon_c] = ratio_and_sd(M.gene.nnon_c,M.gene.Nnon_c);

k = 1.96;
rbagelflank_c = max(0,rbagelflank_c-k*sdbagelflank_c);
rbagelsil_c = max(0,rbagelsil_c-k*sdbagelsil_c);
rflank_c = max(0,rflank_c-k*sdflank_c);
rsil_c = max(0,rsil_c-k*sdsil_c);
rsf_c = max(0,rsf_c-k*sdsf_c);
rsphere_c = max(0,rsphere_c-k*sdsphere_c);
rnon_c = max(0,rnon_c-k*sdnon_c);

rbagelflank_c(M.gene.Nbagelflank_c==0) = 0;
rbagelsil_c(M.gene.Nbagelsil_c==0) = 0;
rflank_c(M.gene.Nflank_c==0) = 0;
rsil_c(M.gene.Nsil_c==0) = 0;
rsf_c(M.gene.Nsf_c==0) = 0;
rsphere_c(M.gene.Nsphere_c==0) = 0;
rnon_c(M.gene.Nnon_c==0) = 0;

fbagelflank_c = bsxfun(@rdivide,rbagelflank_c,M.rflank_c); fbagelflank_c(isnan(fbagelflank_c))=0;
fbagelsil_c = bsxfun(@rdivide,rbagelsil_c,M.rsil_c); fbagelsil_c(isnan(fbagelsil_c))=0;
fflank_c = bsxfun(@rdivide,rflank_c,M.rflank_c); fflank_c(isnan(fflank_c))=0;
fsil_c = bsxfun(@rdivide,rsil_c,M.rsil_c); fsil_c(isnan(fsil_c))=0;
fsf_c = bsxfun(@rdivide,rsf_c,M.rsf_c); fsf_c(isnan(fsf_c))=0;
fsphere_c = bsxfun(@rdivide,rsphere_c,M.rsf_c); fsphere_c(isnan(fsphere_c))=0;
fnon_c = bsxfun(@rdivide,rnon_c,M.rnon_c); fnon_c(isnan(fnon_c))=0;

M.gene.fbagelflank_c = round(fbagelflank_c);
M.gene.fbagelsil_c = round(fbagelsil_c);
M.gene.fflank_c = round(fflank_c);
M.gene.fsil_c = round(fsil_c);
M.gene.fsf_c = round(fsf_c);
M.gene.fsphere_c = round(fsphere_c);
M.gene.fnon_c = round(fnon_c);

%%%%%%%%%%% FINAL SCORING

method = 2;

bkgd = [fbagelflank_c fbagelsil_c fflank_c fsil_c fsf_c fsphere_c];
signal = [fnon_c];

bkgd = max(bkgd,[],2);
signal = max(signal,[],2);

bkgd(bkgd<0.1) = 0.1;
bkgd(isnan(bkgd)) = 1;
bkgd(isinf(bkgd)) = 10;

signal(isnan(signal)) = 0;
signal(isinf(signal)) = 10;

M.gene.bkgd_c = bkgd;
M.gene.signal_c = signal;

score_c = signal./bkgd;

M.gene.score_c = score_c;


%%%%%% SPECIAL CASE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  genes with RPKM=0  --> set score = -100              %
M.gene.score_c(M.V.val{16}==0) = -100;                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% reporting structure
M.g = [];
M.g.name = M.gene.name;

M.g.N = M.gene.Nnon_c(:,end);
M.g.nflank = sum(M.gene.nflank_c,2);
M.g.nsil = sum(M.gene.nsil_c,2);
M.g.nmis = sum(M.gene.nnon_c(:,1:6),2);
M.g.nnull = M.gene.nnon_c(:,7);
M.g.nnull2 = M.gene.nnon_c(:,8);

nflank = sum(M.gene.nflank_c,2); Nflank = M.gene.Nflank_c(:,end);
fflank = (nflank./Nflank)/(sum(nflank)/sum(Nflank)); fflank(Nflank==0) = 0;

nsil = sum(M.gene.nsil_c,2); Nsil = M.gene.Nsil_c(:,end);
fsil = (nsil./Nsil)/(sum(nsil)/sum(Nsil)); fsil(Nsil==0) = 0;

nsf = sum(M.gene.nsf_c,2); Nsf = M.gene.Nsf_c(:,end);
fsf = (nsf./Nsf)/(sum(nsf)/sum(Nsf)); fsf(Nsf==0) = 0;

nbagelflank = sum(M.gene.nbagelflank_c,2); Nbagelflank = M.gene.Nbagelflank_c(:,end);
fbagelflank = (nbagelflank./Nbagelflank)/(sum(nflank)/sum(Nflank)); fbagelflank(Nbagelflank==0) = 0;

nbagelsil = sum(M.gene.nbagelsil_c,2); Nbagelsil = M.gene.Nbagelsil_c(:,end);
fbagelsil = (nbagelsil./Nbagelsil)/(sum(nsil)/sum(Nsil)); fbagelsil(Nbagelsil==0) = 0;

nsf = sum(M.gene.nsf_c,2); Nsf = M.gene.Nsf_c(:,end);
fsf = (nsf./Nsf)/(sum(nsf)/sum(Nsf)); fsf(Nsf==0) = 0;

nsphere = sum(M.gene.nsphere_c,2); Nsphere = M.gene.Nsphere_c(:,end);
fsphere = (nsphere./Nsphere)/(sum(nsf)/sum(Nsf)); fsphere(Nsphere==0) = 0;

nnon = sum(M.gene.nnon_c,2); Nnon = M.gene.Nnon_c(:,end);
fnon = (nnon./Nnon)/(sum(nnon)/sum(Nnon)); fnon(Nnon==0) = 0;

M.g.Fbflk = round(fbagelflank);
M.g.Fbsil = round(fbagelsil);
M.g.Fflk = round(fflank);
M.g.Fsil = round(fsil);
M.g.Fsf = round(fsf);
M.g.Fsphere = round(fsphere);
M.g.Fnon = round(fnon);

M.g.score = M.gene.score_c;
M.g = sort_struct(M.g,'score',-1);





