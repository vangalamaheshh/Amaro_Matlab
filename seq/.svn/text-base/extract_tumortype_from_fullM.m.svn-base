function M = extract_tumortype_from_fullM(F,tti)

if isnan(tti), error('tti is nan'); end

pidx = find(F.pat.ttype_idx==tti);

M=[];
M.ttype = F.ttype.name{tti};
M.gene = F.gene;
M.pat = reorder_struct(F.pat,pidx);
M.cat = F.cat;

ng = slength(M.gene);
np = slength(M.pat);
nc = slength(M.cat);
z = zeros(ng,nc,np);
M.Nsil = z;
M.nsil = z;
M.Nflank = z;
M.nflank = z;
M.Nnon = z;
M.nnon = z;
M.nnon_ignoring_null_categ = z;

cols = (1:nc);  % if double_nulls have been added, this will fortuitously use F's column8=total
gidx = listmap(M.gene.name,F.cov.gene.name);
M.Nsil = round(repmat(F.cov.gene_sil_repcov(gidx,cols),[1 1 np]));       
M.Nnon = round(repmat(F.cov.gene_non_repcov(gidx,cols),[1 1 np]));
M.Nflank = round(repmat(F.cov.gene_flank_repcov(gidx,cols),[1 1 np]));

fprintf('patient:');
for i=1:np, if ~mod(i,10), fprintf(' %d/%d',i,np); end
  midx = find(F.mut.pat_idx==pidx(i));
  midx2 = midx(F.mut.is_silent(midx));
  M.nsil(:,:,i) = hist2d_fast(F.mut.gene_idx(midx2),F.mut.categ(midx2),1,ng,1,nc);
  midx2 = midx(~F.mut.is_silent(midx));
  M.nnon(:,:,i) = hist2d_fast(F.mut.gene_idx(midx2),F.mut.categ(midx2),1,ng,1,nc);
  M.nnon_ignoring_null_categ(:,:,i) = hist2d_fast(F.mut.gene_idx(midx2),F.mut.categ_ignoring_null_categ(midx2),1,ng,1,nc);
  midx = find(F.mut_flank.pat_idx==pidx(i));
  M.nflank(:,:,i) = hist2d_fast(F.mut_flank.gene_idx(midx),F.mut_flank.categ(midx),1,ng,1,nc);
end, fprintf('\n');

M.gene = rmfield_if_exist(M.gene,{'ttype_mutsig_p_ns_s','ttype_mutsig_p','ttype_mutsig_q','bin','binrelrate'});

% remove genes lacking coverage

thresh = 100;  % minimum number of covered bases per gene (total across all patients and categories)
idx = find(sum(sum(M.Nnon,3),2)>thresh & sum(sum(M.Nsil,3),2)>thresh & sum(sum(M.Nflank,3),2)>thresh);
if length(idx)<ng
  fprintf('Removing %d/%d genes that lack coverage\n',ng-length(idx),ng);
  M.Nsil = M.Nsil(idx,:,:);
  M.Nflank = M.Nflank(idx,:,:);
  M.Nnon = M.Nnon(idx,:,:);
  M.nsil = M.nsil(idx,:,:);
  M.nflank = M.nflank(idx,:,:);
  M.nnon = M.nnon(idx,:,:);
  M.nnon_ignoring_null_categ = M.nnon_ignoring_null_categ(idx,:,:);
  M.gene = reorder_struct(M.gene,idx);
end




