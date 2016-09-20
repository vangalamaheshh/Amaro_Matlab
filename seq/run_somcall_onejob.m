function run_somcall_onejob(sample)
% run_somcall_onejob(sample)
%
% ussage example: bsub matlab -nodisplay -r "run_somcall_onejob('gbm/0145/wgs')"
%
t_boomdir = ['/xchip/tcga_scratch/lawrence/' sample '/tumor.boom'];
n_boomdir = ['/xchip/tcga_scratch/lawrence/' sample '/normal.boom'];
outdir = ['/xchip/tcga_scratch/lawrence/' sample '/somcall'];

for chr=1:24
  if exist([outdir '/chr' num2str(chr) '.txt'],'file'), continue; end
  scatter_somcall_chr(t_boomdir,n_boomdir,outdir,chr);
end

% gather

system(['cat ' outdir '/chr* > ' outdir '/exome.txt']);
X = load_struct([outdir '/exome.txt'],'%f%f%s%s%f%f%s%s%s%f',0);
X = rename_fields(X,colx(1:10),{'chr','pos','ref','var','tlod','nlod','fmapqz_TN','avgnmm_TN','nuwp_TN','filter'});
nm = slength(X);

% annotate

M = [];
M.build = repmat({'hg18'},nm,1);
M.chr = X.chr;
M.start = X.pos;
M.end = X.pos;
M.ref = X.ref;
M.tum1 = X.ref;
M.tum2 = X.var;
M.tumbarcode = repmat({[sample '_tumor']},nm,1);
M.normbarcode = repmat({[sample '_normal']},nm,1);

save_struct(M,[outdir '/exome.maf'],'no_headers');
annotate_maflite([outdir '/exome.maf'],[outdir '/exome_annotated.maf']);

M = load_struct([outdir '/exome_annotated.maf']);
M.tumor_LOD = X.tlod;
M.normal_LOD = X.nlod;
M.fmapqz_TN = X.fmapqz_TN;
M.avgnmm_TN = X.avgnmm_TN;
M.nuwp_TN = X.nuwp_TN;
M.filter = X.filter;

save_struct(M,[outdir '/exome_annotated_filtered.maf']);

% stats

nm = slength(M);
nfm = sum(X.filter>=0.5);

S = [];
S = sprintf('%d mutations total\n',nm);
S = [S count_sprintf(M.type) sprintf('\n')];

S = [S sprintf('%d mutations with filter >= 0.5\n',nfm)];
S = [S count_sprintf(M.type(M.filter>=0.5)) sprintf('\n')];

save_textfile(S,[outdir '/exome_stats.txt']);

fprintf('Done.\n');
