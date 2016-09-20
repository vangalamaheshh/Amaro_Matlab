function M = expand_to_n_and_N(M,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'covfile','/xchip/cga/reference/mutsig_params/mel_luad_lusc_mediancov.mat');
% P.covfile is never actually used in this function!


if ~isfield(M.patient,'callscheme')
  if ~isfield(M.patient,'calltype')
    error('Can''t find callscheme/calltype information in M.patient');
  else
    M.patient = rename_field(M.patient,'calltype','callscheme');
  end
end

% mutation tallies
fprintf('Building mutation counts tables: ');
z = nan(M.ng,M.ncat+1,M.np);
M.n_flank = z;
M.n_silent = z;
M.n_nonsilent = z;
M.n_nonsilent_ignoring_null_categ = z;
for c=1:M.ncat, fprintf('%d/%d ',c,M.ncat);
  idx = (M.mut.categ==c);
  idx2 = idx & M.mut.is_flank;
  M.n_flank(:,c,:) = hist2d_fast(M.mut.gene_idx(idx2),M.mut.pat_idx(idx2),1,M.ng,1,M.np);
  idx2 = idx & M.mut.is_coding & M.mut.is_silent;
  M.n_silent(:,c,:) = hist2d_fast(M.mut.gene_idx(idx2),M.mut.pat_idx(idx2),1,M.ng,1,M.np);
  idx2 = idx & M.mut.is_coding & ~M.mut.is_silent;
  M.n_nonsilent(:,c,:) = hist2d_fast(M.mut.gene_idx(idx2),M.mut.pat_idx(idx2),1,M.ng,1,M.np);
  idx = (M.mut.categ_ignoring_null_categ==c);
  idx2 = idx & M.mut.is_coding & ~M.mut.is_silent;
  M.n_nonsilent_ignoring_null_categ(:,c,:) = hist2d_fast(M.mut.gene_idx(idx2),M.mut.pat_idx(idx2),1,M.ng,1,M.np);
end, fprintf('\n');

% coverage tables
fprintf('Building coverage tables: ');
M.N_terr = M.cov.gene_terr;
M.N_flank_terr = M.cov.gene_flank_terr;
M.N_sil_terr = M.cov.gene_sil_terr;
M.N_non_terr = M.cov.gene_non_terr;
M.N_cov = z;
M.N_flank_cov = z;
M.N_sil_cov = z;
M.N_non_cov = z;
classidx = 2-(M.patient.callscheme==2);   % 1=MEL 2=LUNG
for c=1:M.ncat, fprintf('%d/%d ',c,M.ncat);
  M.N_cov(:,c,:) = M.cov.gene_cov(:,classidx,c);
  M.N_flank_cov(:,c,:) = M.cov.gene_flank_cov(:,classidx,c);
  M.N_sil_cov(:,c,:) = M.cov.gene_sil_cov(:,classidx,c);
  M.N_non_cov(:,c,:) = M.cov.gene_non_cov(:,classidx,c);
end, fprintf('\n');
noflank = (M.patient.callscheme==0);
M.N_flank_cov(:,:,noflank)=0;


% totals
fprintf('Computing totals: ');
flds = grep('^(n|N)_',fieldnames(M));
for i=1:length(flds), fprintf('%d/%d ',i,length(flds));
  if flds{i}(1)=='N'
    M.(flds{i})(:,M.TOT,:) = M.(flds{i})(:,M.TOT-1,:);  % copy from indel terr/cov
  else    % n
    M.(flds{i})(:,M.TOT,:) = sum(M.(flds{i})(:,1:M.ncat,:),2);   % sum all categories
  end
end, fprintf('\n');


