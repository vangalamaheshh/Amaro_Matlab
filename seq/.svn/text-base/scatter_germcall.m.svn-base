function scatter_germcall(boomdir,outdir,sample)

% scatter

jobs=[];
for chr=1:24
  if exist([outdir '/chr' num2str(chr) '.txt'],'file'), continue; end
  banner = [sample 'GRMCL' num2str(chr)];
  cmd = ['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "scatter_germcall_chr(''' boomdir ''',''' outdir ''',' num2str(chr) ')"'''''];
  jobs=[jobs;bsub(cmd,banner)];
end

bwait(jobs);

% gather

system(['cat ' outdir '/chr* > ' outdir '/germ.txt']);
X = load_struct([outdir '/germ.txt'],'%f%f%s%s%f%f%s%s%s%f',0);
X = rename_fields(X,colx(1:9),{'chr','pos','ref','var','lod','fmapqz','avgnmm','nuwp','filter'});
nm = slength(X);

keyboard

% annotate

M = [];
M.build = repmat({'hg18'},nm,1);
M.chr = X.chr;
M.start = X.pos;
M.end = X.pos;
M.ref = X.ref;
M.tum1 = X.ref;
M.tum2 = X.var;
M.tumbarcode = repmat({sample},nm,1);
M.normbarcode = repmat({'germline'},nm,1);

save_struct(M,[outdir '/germ.maf'],'no_headers');
annotate_maflite([outdir '/germ.maf'],[outdir '/germ_annotated.maf']);

M = load_struct([outdir '/germ_annotated.maf']);
M.tumor_LOD = X.tlod;
M.normal_LOD = X.nlod;
M.fmapqz_TN = X.fmapqz_TN;
M.avgnmm_TN = X.avgnmm_TN;
M.nuwp_TN = X.nuwp_TN;
M.filter = X.filter;

save_struct(M,[outdir '/germ_annotated_filtered.maf']);

% stats

nm = slength(M);
nfm = sum(X.filter>=0.5);

S = [];
S = sprintf('%d mutations total\n',nm);
S = [S count_sprintf(M.type) sprintf('\n')];

S = [S sprintf('%d mutations with filter >= 0.5\n',nfm)];
S = [S count_sprintf(M.type(M.filter>=0.5)) sprintf('\n')];

save_textfile(S,[outdir '/germ_stats.txt']);

fprintf('Done.\n');
