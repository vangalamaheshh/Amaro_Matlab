function scatter_somcall(t_boomdir,n_boomdir,outdir,sample,params)
% scatter_somcall(t_boomdir,n_boomdir,outdir,sample,params)
%
% call somatic mutations in exome
%
% Mike Lawrence 2009

if ~exist('sample','var'), sample = ''; end   % used only for LSF jobnames

if ~exist('params','var'), params=[]; end
params=impose_default_value(params,'lane_blacklist_file',[]);
params=impose_default_value(params,'lanetable_file',[]);

if ~exist(t_boomdir,'file'), error('Not found: %s',t_boomdir); end
if ~exist(n_boomdir,'file'), error('Not found: %s',n_boomdir); end
if ~exist(outdir,'dir'), mkdir(outdir); end

if ~isempty(params.lane_blacklist_file)
  fprintf('\n*** WARNING ***  blacklist usage is untested!\n\n');
  subfprintf('Checking reads against lane blacklist\n');
  if isempty(params.lanetable_file), error('Must supply lanetable_file'); end
  LT = load_struct(params.lanetable_file);
  fcl = regexprep(LT.PU,'^(.....)[^\.]*(\.\d)$','$1$2');
  BL = load_lines(params.lane_blacklist_file);
  params.readgroup_blacklist = find(ismember(fcl,BL))-1;
end

params_file=[ outdir '/params.mat'];
save(params_file,'params');

% scatter

jobs=[];
for chr=1:24
  if exist([outdir '/chr' num2str(chr) '.txt'],'file'), continue; end
  banner = [sample 'SOMCL' num2str(chr)];
  cmd = ['-R "rusage[matlab=1:mem=8]" ''''matlab -nodisplay '...
         '-r "scatter_somcall_chr(''' t_boomdir ''',''' n_boomdir ''',''' outdir ''',' num2str(chr) ...
         ',''' params_file ...
         ''')"'''''];
  disp(cmd);
  jobs=[jobs;bsub(cmd,banner)];
end

bwait(jobs);

% gather

M = cell(24,1);
for chr=1:24
  tmp = load([outdir '/chr' num2str(chr) '.mat']);
  M{chr} = tmp.M;
end
M = cat(1,M{:});   % works OK with emptys

% annotate

T = annotate_M_calls(M);
save_struct(T,[outdir '/exome.txt']);

% stats

T = make_numeric(T,{'filter'});
nm = slength(T);
nfm = sum(T.filter>=0.5);

S = [];
S = sprintf('%d mutations total\n',nm);
S = [S count_sprintf(T.type) sprintf('\n')];

S = [S sprintf('%d mutations with filter >= 0.5\n',nfm)];
S = [S count_sprintf(T.type(T.filter>=0.5)) sprintf('\n')];

save_textfile(S,[outdir '/exome_stats.txt']);

fprintf('Done.\n');

